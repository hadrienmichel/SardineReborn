"""Sardine Reborn desktop application for seismic-refraction processing.

The module provides the complete Qt user interface and the numerical helpers
used to load SEG-2 data, pick first arrivals, build layered starting models,
run pyGIMLi travel-time inversions, and display or export the results.

Times are expressed in seconds, distances in metres, and velocities in metres
per second unless a function explicitly states otherwise.
"""

# TODO/IDEAS
# Add option to merge similar sources/receivers at loading of geometry file (new files to add with similar paths (roll-along support))
# Add posibilities for auto-picking (DONE on 04-09-2026 David Caterina)
# Debug Modelling with multiple datasets (To Test)
# Add option to visualize the FFT of the datasets (DONE on 16-05-2023)
# Add option to pick along a line (DONE on 15-05-2023)
# Add option to change the header (dt for example)
# Add options for setting t0 : by value, by graphical picking, by line picking (automated?) (DONE on 22-05-2023) (ADDED an option to set t0 for the entire data based on geometry 04-09-2026 David Caterina)
# If unable to pick negative times --> impossible to change offset when issue... (DONE on 04-09-2026 David Caterina)
# If offset from station --> issue with graph in set-t0
# The offset datasets are currently still openning (bug) the set t=0 window with the graphs, even though the correct behaviour would be to ask to pick some traces befor hand.

# Imports for the inner functions
import sys
import os
import re
from copy import deepcopy
import typing
import numpy as np
# pyGIMLi <= 1.3 still uses NumPy's deprecated scalar aliases. They were
# removed in NumPy 1.24, so restore only the aliases required by that legacy
# release before importing pyGIMLi. Newer pyGIMLi versions ignore this shim.
_NUMPY_LEGACY_ALIASES = {
    'float': float,
    'int': int,
    'bool': bool,
    'complex': complex,
}
for _alias, _builtin_type in _NUMPY_LEGACY_ALIASES.items():
    if _alias not in np.__dict__:
        setattr(np, _alias, _builtin_type)
import time
import matplotlib
from scipy import stats
from scipy.signal import find_peaks
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib import animation
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.backend_bases import MouseButton
from matplotlib import pyplot
from matplotlib import font_manager
# Imports for the seismic data input
from obspy import read
try:
    from obspy.signal.trigger import aic_simple, classic_sta_lta
    AIC_SIMPLE_IMPORT = True
except ImportError:
    AIC_SIMPLE_IMPORT = False
# Imports for the data inversion
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics import TravelTimeManager as TTMgr
from pygimli.physics.traveltime import drawFirstPicks
try:
    from pygimli.physics.traveltime.ratools import createGradientModel2D
except ImportError:
    from pygimli.physics.traveltime import createGradientModel2D
from pygimli.viewer import show as pgshow
# Imports for saving/reading states:
import pickle
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import (
    QApplication,
    QWidget,
    QTabWidget,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QAction,
    QMainWindow,
    QStatusBar,
    QMessageBox,
    QSlider,
    QDialog,
    QFileDialog,
    QInputDialog,
    QProgressBar,
    QComboBox,
    QPushButton,
    QCheckBox,
    QSpinBox,
    QFontComboBox,
    QDoubleSpinBox,
    QLabel,
    QGroupBox,
    QLineEdit,
    QShortcut)
from PyQt5 import QtCore, QtGui
from PyQt5.QtGui import QKeySequence


matplotlib.use('Qt5Agg')
NavigationToolbar2QT.toolitems = [('Home', 'Reset original view', 'home', 'home'),
                                  (None, None, None, None),
                                  ('Pan', 'Left button pans, Ri...xes aspect',
                                   'move', 'pan'),
                                  ('Zoom', 'Zoom to rectangle\nx/...xes aspect',
                                   'zoom_to_rect', 'zoom'),
                                  (None, None, None, None),
                                  ('Save', 'Save the figure', 'filesave', 'save_figure')]

DEFAULT_STATUS = "Idle."
DEFAULT_ERROR = 0.03  # By default, the error on the picking is going to be 3%


def absolute_pick_error(relative_error, travel_time, sample_interval):
    """Return a strictly positive travel-time uncertainty in seconds.

    Parameters
    ----------
    relative_error : float
        Relative picking uncertainty expressed as a fraction of travel time.
        Invalid or non-positive values are replaced by ``DEFAULT_ERROR``.
    travel_time : float
        Picked travel time in seconds.
    sample_interval : float
        Sampling interval in seconds.

    Returns
    -------
    float
        Absolute uncertainty in seconds, with a floor of half a sample or one
        microsecond, whichever is larger.

    Notes
    -----
    The picker stores relative errors internally.  pyGIMLi's ``err`` field,
    however, is an absolute uncertainty in the same unit as ``t``.  A
    half-sample floor keeps the exported uncertainty meaningful at (or before)
    t=0 and protects the inversion from null or negative values.
    """
    try:
        relative_error = float(relative_error)
    except (TypeError, ValueError):
        relative_error = DEFAULT_ERROR
    if not np.isfinite(relative_error) or relative_error <= 0:
        relative_error = DEFAULT_ERROR
    try:
        travel_time = float(travel_time)
    except (TypeError, ValueError):
        travel_time = 0.0
    try:
        sample_interval = abs(float(sample_interval))
    except (TypeError, ValueError):
        sample_interval = 0.0
    minimum_error = max(0.5*sample_interval, 1e-6)
    return max(abs(relative_error)*abs(travel_time), minimum_error)


def build_model(sourceX, receiversX, times, nbLayers=2, orientation=1):
    """Estimate a layered starting model from one travel-time branch.

    Parameters
    ----------
    sourceX : float
        Horizontal source coordinate in metres.
    receiversX : array-like
        Horizontal receiver coordinates in metres.
    times : array-like
        First-arrival travel times in seconds, ordered like ``receiversX``.
    nbLayers : int, default=2
        Number of linear hodograph segments to estimate.
    orientation : {1, -1}, default=1
        Profile direction relative to the source.

    Returns
    -------
    intercept_times : numpy.ndarray
        Time-axis intercept of every fitted hodograph segment, in seconds.
    apparent_velocities : numpy.ndarray
        Apparent velocity of every segment, in metres per second.
    points : numpy.ndarray
        ``(x, time)`` coordinates delimiting the fitted segments.

    Notes
    -----
    Curvature maxima define candidate changes of slope. The result is intended
    only as an inversion starting model, not as a final geological model.
    """
    # Change the coordinates if requiered and normalize to sourceX = 0:
    receiversXSave = receiversX
    if orientation < 0:
        receiversX = np.flip(np.abs(receiversX - sourceX))
        times = np.flip(times)
    else:
        receiversX = receiversX - sourceX
    curvature = np.abs(np.gradient(np.gradient(times, receiversX), receiversX))
    # Find the highest curvatures in the set:
    maxCurvIndex = np.sort(np.argpartition(
        curvature, -(nbLayers-1))[-(nbLayers-1):])
    maxCurvIndex = np.concatenate(([0], maxCurvIndex, [len(times)-1]))
    v = np.zeros((nbLayers,))
    inter = np.zeros((nbLayers,))
    # Points to store the intersections of interest
    points = np.zeros((nbLayers+1, 2))
    for i in range(nbLayers):
        x = receiversX[maxCurvIndex[i]: maxCurvIndex[i+1]]
        y = times[maxCurvIndex[i]: maxCurvIndex[i+1]]
        if i == 0:
            x = x[:, np.newaxis]
            p = np.linalg.lstsq(x, y, rcond=None)
            v[i] = 1/p[0][0]
            inter[i] = 0
        else:
            p = np.polyfit(x, y, 1)
            v[i] = 1/p[0]
            inter[i] = p[1]

    points[0, 0] = sourceX
    for i in range(nbLayers):
        if i < nbLayers-1:
            xTemp = (inter[i]-inter[i+1])/((1/v[i+1])-(1/v[i]))
            tTemp = inter[i] + xTemp*(1/v[i])
            points[i+1, 0] = sourceX + orientation*xTemp
            points[i+1, 1] = tTemp
        else:
            if orientation > 0:
                points[i+1, 0] = max(receiversXSave)
            else:
                points[i+1, 0] = min(receiversXSave)
            points[i+1, 1] = inter[-1] + \
                np.abs(sourceX - points[i+1, 0])*(1/v[-1])

    return inter, v, points


def model1D(inter, v):
    """Derive horizontal-layer thicknesses from intercept times and velocities.

    Parameters
    ----------
    inter : array-like
        Intercept time for each refracted branch, in seconds.
    v : array-like
        Layer velocities in metres per second, ordered from top to bottom.

    Returns
    -------
    thicknesses : numpy.ndarray
        Estimated thickness of each layer above the half-space, in metres.
    velocities : array-like
        The input velocity vector, returned unchanged.
    positions : numpy.ndarray
        Cumulative horizontal and depth coordinates of the interfaces.
    """
    h = np.ones((len(v)-1,))
    pos = np.zeros((len(v), 2))
    for i in np.arange(1, len(v)):
        iCr = np.arcsin(v[:i]/v[i])
        tInterp = inter[i]
        for j in np.arange(i-1):
            tInterp -= 2*h[j]*np.cos(iCr[j])/v[j]
            if j < i-1:
                pos[i, 0] += np.tan(iCr[j]) * h[j]
        h[i-1] = v[i-1] * tInterp / (2*np.cos(iCr[i-1]))
        pos[i, 1] = pos[i-1, 1] + h[i-1]
        pos[i, 0] = pos[i-1, 0] + np.tan(iCr[i-1]) * h[i-1]
    return h, v, pos


def modelWithSlope(interS, vS):
    """Estimate dipping-layer geometry from reciprocal apparent velocities.

    Parameters
    ----------
    interS : array-like
        Intercept times for forward and reverse shots, with one row per layer.
    vS : array-like
        Forward and reverse apparent velocities in metres per second, with one
        row per layer.

    Returns
    -------
    velocities : list of float
        Estimated true layer velocities in metres per second.
    left_thicknesses, right_thicknesses : list of float
        Interface-normal thicknesses below the two profile ends, in metres.

    Notes
    -----
    The calculation follows Mota (1954), "Determination of dips and depths of
    geological layers by the seismic refraction method", for two or three
    layers.
    """
    nbLayers = np.shape(vS)[0]
    # We take the mean velocity between the two possibilities
    v0 = np.sum(vS[0, :])/2
    # Find alpha and v2
    ipA = np.arcsin(v0/vS[1, 0])
    imA = np.arcsin(v0/vS[1, 1])
    iCr = (ipA + imA)/2
    theta1 = (imA - ipA)/2
    v1 = v0/np.sin(iCr)
    zL = [v0*interS[1, 0]/(2*np.cos(iCr))]
    zR = [v0*interS[1, 1]/(2*np.cos(iCr))]
    hL = [zL[0]/np.cos(theta1)]  # Depth of the interface below the source A
    hR = [zR[0]/np.cos(theta1)]  # Depth of the interface below the source B
    v = [v0, v1]
    if nbLayers > 2:
        alpha = np.arcsin(v0/vS[2, 1]) - theta1
        beta = np.arcsin(v0/vS[2, 0]) + theta1
        gamma = np.arcsin(v1/v0 * np.sin(alpha))
        delta = np.arcsin(v1/v0 * np.sin(beta))
        iCr = (gamma+delta)/2
        theta2 = (gamma-delta)/2 + theta1
        v2 = v1/np.sin(iCr)
        zL.append(v1*(interS[2, 0] - zL[0]/v0 *
                  (np.cos(alpha+beta)+1)/np.cos(alpha))/(2*np.cos(iCr)))
        zR.append(v1*(interS[2, 1] - zR[0]/v0 *
                  (np.cos(alpha+beta)+1)/np.cos(beta))/(2*np.cos(iCr)))
        hL.append(1/np.cos(theta2) *
                  (zL[0]*np.cos(alpha-theta2+theta1)/np.cos(alpha) + zL[1]))
        hR.append(1/np.cos(theta2) *
                  (zR[0]*np.cos(beta+theta2-theta1)/np.cos(beta) + zR[1]))
        v = [v0, v1, v2]
    return v, hL, hR


def calculateDistance(pts, pt):
    """Return the Euclidean distance from every point in ``pts`` to ``pt``.

    Parameters
    ----------
    pts : array-like, shape (n, 2)
        Coordinates of the points whose distances are required.
    pt : array-like, shape (2,)
        Reference coordinate.

    Returns
    -------
    numpy.ndarray, shape (n,)
        Distance from each point to the reference coordinate.
    """
    pts = np.asarray(pts)
    xDiff = pts[:, 0] - pt[0]
    yDiff = pts[:, 1] - pt[1]
    dist = np.sqrt(xDiff**2 + yDiff**2)
    return dist


def hybrid_pick_candidates(signal, start_idx, stop_idx, max_candidates=6):
    """Return AIC-refined STA/LTA candidates for one seismic trace.

    Parameters
    ----------
    signal : array-like
        One-dimensional seismic-amplitude samples.
    start_idx, stop_idx : int
        Half-open sample-index interval in which candidates are sought.
    max_candidates : int, default=6
        Maximum number of candidates returned after duplicate removal.

    Returns
    -------
    list of dict
        Candidates sorted by decreasing confidence. Each dictionary contains
        ``index``, ``confidence``, ``cost``, and ``snr``. An empty list means
        that no reliable onset could be extracted.
    """
    values = np.asarray(signal, dtype=float)
    start_idx = max(0, int(start_idx))
    stop_idx = min(values.size, int(stop_idx))
    segment = values[start_idx:stop_idx]
    if segment.size < 24:
        return []

    finite = np.isfinite(segment)
    if not finite.any():
        return []
    fill_value = float(np.nanmedian(segment[finite]))
    segment = np.where(finite, segment, fill_value)
    segment = segment - np.median(segment)
    mad = np.median(np.abs(segment)) / 0.67448975
    scale = mad if np.isfinite(mad) and mad > 0 else np.std(segment)
    if not np.isfinite(scale) or scale <= np.finfo(float).eps:
        return []
    normalized = np.clip(segment / scale, -25.0, 25.0)

    npts = normalized.size
    nsta = max(3, min(24, npts // 80))
    nlta = max(5 * nsta, min(240, npts // 8))
    if nlta >= npts - 4:
        nlta = max(nsta + 2, npts // 3)
    characteristic = np.asarray(
        classic_sta_lta(normalized, nsta, nlta), dtype=float)
    characteristic[~np.isfinite(characteristic)] = 0.0

    baseline = characteristic[nlta:] if nlta < npts else characteristic
    baseline_median = float(np.median(baseline)) if baseline.size else 0.0
    baseline_mad = float(np.median(np.abs(baseline-baseline_median))) if baseline.size else 0.0
    threshold = max(2.0, baseline_median + 5.0 * baseline_mad)

    above = characteristic >= threshold
    crossings = np.flatnonzero(above & ~np.r_[False, above[:-1]])
    peaks, _ = find_peaks(characteristic, distance=max(2, 2*nsta))
    if peaks.size:
        strongest = peaks[np.argsort(characteristic[peaks])[-max_candidates:]]
    else:
        strongest = np.array([], dtype=int)

    global_aic = np.asarray(aic_simple(normalized), dtype=float)
    global_aic[~np.isfinite(global_aic)] = np.inf
    aic_minima, _ = find_peaks(-global_aic, distance=max(3, 3*nsta))
    if aic_minima.size:
        best_aic = aic_minima[np.argsort(global_aic[aic_minima])[:max_candidates]]
    else:
        best_aic = np.array([int(np.argmin(global_aic))])

    seeds = np.unique(np.r_[crossings[:max_candidates], strongest, best_aic])
    candidates = []
    radius = max(20, 5*nsta)
    for seed in seeds:
        local_start = max(0, int(seed)-radius)
        local_stop = min(npts, int(seed)+radius+1)
        local = normalized[local_start:local_stop]
        if local.size < 8:
            continue
        local_aic = np.asarray(aic_simple(local), dtype=float)
        margin = max(2, int(round(0.05*local_aic.size)))
        valid_aic = local_aic[margin:local_aic.size-margin]
        if valid_aic.size == 0 or not np.isfinite(valid_aic).any():
            continue
        refined = local_start + margin + int(np.nanargmin(valid_aic))

        pre = normalized[max(0, refined-nlta):refined]
        post = normalized[refined:min(npts, refined+3*nsta)]
        pre_rms = np.sqrt(np.mean(pre**2)) if pre.size else 1.0
        post_rms = np.sqrt(np.mean(post**2)) if post.size else 0.0
        snr = post_rms / max(pre_rms, np.finfo(float).eps)
        cf_start = max(0, refined-nsta)
        cf_stop = min(npts, refined+2*nsta+1)
        strength = float(np.max(characteristic[cf_start:cf_stop]))
        finite_aic = valid_aic[np.isfinite(valid_aic)]
        aic_scale = np.std(finite_aic)
        aic_contrast = ((np.median(finite_aic)-np.min(finite_aic)) /
                        max(aic_scale, np.finfo(float).eps))

        strength_score = np.clip(
            (strength-1.0)/max(threshold-1.0, 1.0), 0.0, 1.0)
        # Preserve discrimination between a weak secondary change and the very
        # sharp energy increase usually associated with the first break.
        snr_score = np.clip(np.log(max(snr, 1.0))/np.log(50.0), 0.0, 1.0)
        aic_score = np.clip(aic_contrast/4.0, 0.0, 1.0)
        confidence = float(
            0.45*strength_score + 0.35*snr_score + 0.20*aic_score)
        if not np.isfinite(confidence):
            confidence = 0.0
        candidates.append({
            'index': start_idx + refined,
            'confidence': confidence,
            'cost': 1.0-confidence,
            'snr': float(snr)})

    # Merge candidates that converged to the same local AIC minimum.
    unique = {}
    for candidate in candidates:
        idx = candidate['index']
        if idx not in unique or candidate['confidence'] > unique[idx]['confidence']:
            unique[idx] = candidate
    candidates = sorted(unique.values(), key=lambda item: item['confidence'], reverse=True)
    return candidates[:max_candidates]


def near_source_energy_candidate(signal, begin_time, dt,
                                 start_idx, stop_idx):
    """Detect the first persistent energy rise close to a seismic source.

    Parameters
    ----------
    signal : array-like
        One-dimensional seismic-amplitude samples.
    begin_time : float
        Time of the first sample in seconds, relative to the current t0.
    dt : float
        Sampling interval in seconds.
    start_idx, stop_idx : int
        Half-open sample-index interval allowed for the search.

    Returns
    -------
    dict or None
        Candidate metadata containing its sample ``index``, ``confidence``,
        ``cost``, ``snr``, and a ``near_source_onset`` flag; ``None`` when no
        persistent high-SNR energy onset satisfies the checks.
    """
    values = np.asarray(signal, dtype=float)
    start_idx = max(0, int(start_idx))
    stop_idx = min(values.size, int(stop_idx))
    if stop_idx-start_idx < 64:
        return None
    finite = np.isfinite(values)
    if not finite.any():
        return None
    values = np.where(finite, values, np.nanmedian(values[finite]))

    zero_idx = int(round((0.0-begin_time)/dt))
    baseline_start = max(start_idx, zero_idx-int(round(0.030/dt)))
    baseline_stop = min(stop_idx, zero_idx-int(round(0.005/dt)))
    if baseline_stop-baseline_start < 32:
        baseline_stop = min(
            stop_idx, start_idx+max(32, (stop_idx-start_idx)//6))
    if baseline_stop-baseline_start < 16:
        return None
    values = values-np.median(values[baseline_start:baseline_stop])

    smooth_samples = max(3, int(round(0.0005/dt)))
    energy = np.sqrt(np.convolve(
        values*values, np.ones(smooth_samples)/smooth_samples,
        mode='same'))
    baseline = energy[baseline_start:baseline_stop]
    centre = float(np.median(baseline))
    spread = 1.4826*float(np.median(np.abs(baseline-centre)))
    if not np.isfinite(spread) or spread <= np.finfo(float).eps:
        spread = float(np.std(baseline))
    if not np.isfinite(spread) or spread <= np.finfo(float).eps:
        return None

    search_start = max(start_idx, zero_idx-int(round(0.015/dt)))
    search_stop = min(stop_idx, zero_idx+int(round(0.030/dt)))
    if search_stop <= search_start:
        return None
    above = energy[search_start:search_stop] >= centre+6.0*spread
    persistence = max(2, int(round(0.00025/dt)))
    if above.size < persistence:
        return None
    persistent = np.convolve(
        above.astype(int), np.ones(persistence, dtype=int), mode='valid')
    starts = np.flatnonzero(persistent >= persistence)
    if not starts.size:
        return None
    onset = search_start+int(starts[0])

    post_stop = min(values.size, onset+max(
        smooth_samples, int(round(0.002/dt))))
    pre_start = max(start_idx, onset-int(round(0.005/dt)))
    pre = values[pre_start:onset]
    post = values[onset:post_stop]
    pre_rms = np.sqrt(np.mean(pre**2)) if pre.size else centre
    post_rms = np.sqrt(np.mean(post**2)) if post.size else 0.0
    snr = post_rms/max(pre_rms, np.finfo(float).eps)
    if not np.isfinite(snr) or snr < 3.0:
        return None
    peak_strength = ((np.max(energy[onset:post_stop])-centre) /
                     max(spread, np.finfo(float).eps))
    confidence = float(np.clip(
        0.68+0.12*np.log10(max(peak_strength, 1.0)), 0.68, 0.97))
    return {
        'index': onset,
        'confidence': confidence,
        'cost': 1.0-confidence,
        'snr': float(snr),
        'near_source_onset': True}


def _robust_line_prediction(x, y, query_x=None):
    """Fit a line without letting one late phase control the prediction.

    Parameters
    ----------
    x, y : array-like
        Coordinates of the observations used to fit the trend.
    query_x : array-like or None, optional
        Coordinates at which to evaluate the trend. The values in ``x`` are
        used when this argument is omitted.

    Returns
    -------
    numpy.ndarray
        Predicted ``y`` values at ``query_x``. A Theil-Sen fit is used for at
        least three observations; smaller inputs use a linear fit.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    query_x = x if query_x is None else np.asarray(query_x, dtype=float)
    if y.size == 1 or np.ptp(x) <= np.finfo(float).eps:
        return np.full(query_x.shape, float(np.median(y)))
    if y.size >= 3:
        estimate = stats.theilslopes(y, x)
        slope = estimate.slope if hasattr(estimate, 'slope') else estimate[0]
        intercept = (estimate.intercept if hasattr(estimate, 'intercept')
                     else estimate[1])
    else:
        slope, intercept = np.polyfit(x, y, 1)
    return intercept+slope*query_x


def _candidate_for_trend(candidates, prediction, tolerance):
    """Return the best waveform candidate close enough to one trend.

    Parameters
    ----------
    candidates : sequence of dict
        Waveform candidates containing at least ``index`` and ``cost``.
    prediction : float
        Expected arrival sample derived from the spatial trend.
    tolerance : float
        Nominal admissible residual in samples.

    Returns
    -------
    tuple
        ``(candidate, score)`` for the lowest-cost nearby candidate, or
        ``(None, numpy.inf)`` if none lies within three tolerances.
    """
    nearby = [candidate for candidate in candidates
              if abs(candidate['index']-prediction) <= 3.0*tolerance]
    if not nearby:
        return None, np.inf
    candidate = min(
        nearby,
        key=lambda item: item['cost'] + 0.03*abs(
            item['index']-prediction)/tolerance)
    score = (candidate['cost'] +
             0.03*abs(candidate['index']-prediction)/tolerance)
    return candidate, float(score)


def _forward_trend_guard(selected, candidates_by_trace, order, x,
                         tolerance):
    """Keep an isolated bad trace out of the trend used by later receivers.

    Parameters
    ----------
    selected : dict
        Current mapping from trace identifiers to selected candidates.
    candidates_by_trace : sequence or mapping
        Candidate lists indexed by trace identifier.
    order : sequence of int
        Trace identifiers ordered away from the source.
    x : array-like
        Spatial coordinate corresponding to each entry in ``order``.
    tolerance : float
        Maximum nominal trend residual in samples.

    Returns
    -------
    tuple
        Updated selections, per-trace forward predictions, and the set of
        traces still suspected of breaking the current trend.
    """
    selected = dict(selected)
    reliable_positions = []
    pending_positions = []
    guarded_predictions = {}
    suspected_traces = set()
    for position, trace_id in enumerate(order):
        candidates = candidates_by_trace[trace_id]
        locked = (len(candidates) == 1 and
                  candidates[0].get('cost', 0.0) < 0)
        if locked or len(reliable_positions) < 3:
            reliable_positions.append(position)
            guarded_predictions[trace_id] = float(selected[trace_id]['index'])
            pending_positions = []
            continue

        history = reliable_positions[-6:]
        history_y = np.asarray(
            [selected[order[i]]['index'] for i in history], dtype=float)
        predicted = float(_robust_line_prediction(
            x[history], history_y, np.asarray([x[position]]))[0])
        replacement, _ = _candidate_for_trend(
            candidates, predicted, tolerance)
        if replacement is not None:
            selected[trace_id] = replacement
            guarded_predictions[trace_id] = predicted
            reliable_positions.append(position)
            pending_positions = []
            continue

        # Do not let this trace update the extrapolation. Three consecutive
        # misses may, however, describe a real new branch and are admitted if
        # they form a coherent line of their own.
        pending_positions.append(position)
        guarded_predictions[trace_id] = predicted
        suspected_traces.add(trace_id)
        if len(pending_positions) >= 3:
            pending = pending_positions[-3:]
            pending_y = np.asarray(
                [selected[order[i]]['index'] for i in pending], dtype=float)
            pending_prediction = _robust_line_prediction(x[pending], pending_y)
            coherent = np.max(np.abs(
                pending_y-pending_prediction)) <= 1.5*tolerance
            confident = np.median([
                selected[order[i]].get('confidence', 0.0)
                for i in pending]) >= 0.40
            if coherent and confident:
                reliable_positions.extend(pending)
                for i, value in zip(pending, pending_prediction):
                    guarded_predictions[order[i]] = float(value)
                    suspected_traces.discard(order[i])
                pending_positions = []
    return selected, guarded_predictions, suspected_traces


def _piecewise_linear_prediction(x, y, base_tolerance=8.0,
                                 min_segment_length=4,
                                 breakpoint_penalty=4.5):
    """Fit a robust piecewise-linear trend to ordered picks.

    Parameters
    ----------
    x, y : array-like
        Ordered spatial coordinates and picked sample indices.
    base_tolerance : float, default=8.0
        Residual scale in samples used by the robust loss.
    min_segment_length : int, default=4
        Minimum number of picks required on each fitted segment.
    breakpoint_penalty : float, default=4.5
        Cost added for every change of slope.

    Returns
    -------
    prediction : numpy.ndarray
        Fitted value corresponding to every input observation.
    breakpoints : list of int
        Indices at which a new linear segment begins.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    count = y.size
    if count == 0:
        return np.array([], dtype=float), []
    if count < 2*min_segment_length:
        return _robust_line_prediction(x, y), []

    scale = max(float(base_tolerance), np.finfo(float).eps)
    segment_cache = {}

    def segment(start, stop):
        """Return the robust fit and normalized cost for one index interval."""
        key = (start, stop)
        if key not in segment_cache:
            prediction = _robust_line_prediction(x[start:stop], y[start:stop])
            residual = np.abs(y[start:stop]-prediction)/scale
            # Huber-like loss: large isolated errors remain bounded enough not
            # to manufacture a false change of slope.
            loss = np.where(residual <= 1.5,
                            0.5*residual**2,
                            1.5*residual-1.125)
            segment_cache[key] = (float(np.sum(loss)), prediction)
        return segment_cache[key]

    best = np.full(count+1, np.inf)
    previous = np.full(count+1, -1, dtype=int)
    best[0] = -breakpoint_penalty  # The first segment is not a breakpoint.
    for stop in range(min_segment_length, count+1):
        starts = [0]
        starts.extend(range(min_segment_length,
                            stop-min_segment_length+1))
        for start in starts:
            if start and not np.isfinite(best[start]):
                continue
            cost = best[start]+segment(start, stop)[0]+breakpoint_penalty
            if cost < best[stop]:
                best[stop] = cost
                previous[stop] = start

    if previous[count] < 0:
        return _robust_line_prediction(x, y), []
    intervals = []
    stop = count
    while stop > 0:
        start = int(previous[stop])
        intervals.append((start, stop))
        stop = start
    intervals.reverse()

    prediction = np.empty(count, dtype=float)
    for start, stop in intervals:
        prediction[start:stop] = segment(start, stop)[1]
    breakpoints = [stop for _, stop in intervals[:-1]]
    return prediction, breakpoints


def select_consistent_candidates(candidates_by_trace, trace_order,
                                 positions=None, base_tolerance=8.0,
                                 regularize_piecewise=True):
    """Select first-arrival candidates with robust spatial coherence.

    Parameters
    ----------
    candidates_by_trace : sequence or mapping
        Candidate dictionaries grouped by trace identifier.
    trace_order : sequence of int
        Trace identifiers ordered away from the source along one branch.
    positions : array-like or None, optional
        Receiver coordinates indexed by trace identifier. Sequential positions
        are used when omitted.
    base_tolerance : float, default=8.0
        Nominal admissible residual between a candidate and trend, in samples.
    regularize_piecewise : bool, default=True
        If true, validate changes of slope with a segmented spatial trend.

    Returns
    -------
    dict
        Mapping from accepted trace identifiers to selected candidate
        dictionaries. Returned candidates also contain spatial-prediction and
        segment metadata when piecewise regularization is enabled.
    """
    independent = {}
    order = []
    for item in trace_order:
        trace_id = int(item)
        candidates = candidates_by_trace[trace_id]
        if not candidates:
            continue
        if (len(candidates) == 1 and
                candidates[0].get('exclude_from_spatial_trend', False)):
            independent[trace_id] = dict(candidates[0])
        else:
            order.append(trace_id)
    if not order:
        return independent

    raw = [min(candidates_by_trace[i], key=lambda item: item['cost'])['index']
           for i in order]
    differences = np.diff(raw)
    expected_step = float(np.median(differences)) if differences.size else 0.0
    deviation = (np.median(np.abs(differences-expected_step))
                 if differences.size else 1.0)
    transition_scale = max(4.0, abs(expected_step)*0.75, 3.0*deviation)

    accumulated = [np.asarray(
        [candidate['cost'] for candidate in candidates_by_trace[order[0]]])]
    backtrack = []
    for position in range(1, len(order)):
        previous = candidates_by_trace[order[position-1]]
        current = candidates_by_trace[order[position]]
        current_cost = np.full(len(current), np.inf)
        current_back = np.zeros(len(current), dtype=int)
        for j, candidate in enumerate(current):
            transitions = []
            for k, previous_candidate in enumerate(previous):
                step = candidate['index']-previous_candidate['index']
                smoothness = 0.30*abs(step-expected_step)/transition_scale
                # Use a bounded robust penalty: a locally bad trace must not
                # drag every subsequent pick towards the same late phase.
                smoothness = min(smoothness, 0.75)
                transitions.append(accumulated[-1][k] + smoothness)
            current_back[j] = int(np.argmin(transitions))
            current_cost[j] = candidate['cost'] + transitions[current_back[j]]
        accumulated.append(current_cost)
        backtrack.append(current_back)

    selected_ids = [int(np.argmin(accumulated[-1]))]
    for links in reversed(backtrack):
        selected_ids.append(int(links[selected_ids[-1]]))
    selected_ids.reverse()
    selected = {trace_id: candidates_by_trace[trace_id][candidate_id]
                for trace_id, candidate_id in zip(order, selected_ids)}

    # Close to the source, a sharp persistent energy onset is more reliable
    # than extrapolating a far-offset phase back through the near field. Keep
    # these signal-validated candidates as the starting geometry for the
    # branch. Manual anchors are unaffected because they do not carry this tag.
    for trace_id in order:
        near_source = [candidate
                       for candidate in candidates_by_trace[trace_id]
                       if candidate.get('near_source_onset', False) and
                       candidate.get('snr', 0.0) >= 3.0]
        if near_source:
            selected[trace_id] = max(
                near_source, key=lambda item: item['confidence'])
    if not regularize_piecewise:
        selected = {trace_id: dict(candidate)
                    for trace_id, candidate in selected.items()}
        selected.update(independent)
        return selected

    # Reconsider every trace against a global segmented trend. A breakpoint is
    # only introduced when at least four consecutive receivers support each
    # side, which implements a short look-ahead and prevents one bad pick from
    # changing the trajectory followed by all subsequent traces.
    if positions is None:
        x = np.arange(len(order), dtype=float)
    else:
        position_values = np.asarray(positions, dtype=float)
        x = np.asarray([position_values[i] for i in order], dtype=float)
    min_segment = 4
    selected, forward_predictions, suspected_traces = _forward_trend_guard(
        selected, candidates_by_trace, order, x, base_tolerance)
    for _ in range(2):
        y = np.asarray([selected[i]['index'] for i in order], dtype=float)
        prediction, breakpoints = _piecewise_linear_prediction(
            x, y, base_tolerance=base_tolerance,
            min_segment_length=min_segment)
        changed = False
        for position, trace_id in enumerate(order):
            candidates = candidates_by_trace[trace_id]
            if len(candidates) == 1 and candidates[0].get('cost', 0.0) < 0:
                continue  # Explicit manual/validated anchor.
            predicted = prediction[position]
            replacement = min(
                candidates,
                key=lambda item: item['cost'] + 0.55*min(
                    abs(item['index']-predicted)/base_tolerance, 2.5))
            forward_prediction = forward_predictions.get(trace_id)
            if (forward_prediction is not None and
                    abs(replacement['index']-forward_prediction) >
                    2.5*base_tolerance):
                continue
            if replacement['index'] != selected[trace_id]['index']:
                selected[trace_id] = replacement
                changed = True
        if not changed:
            break

    # A coherent secondary phase can otherwise manufacture a false breakpoint:
    # one unusable trace starts it and the following traces are pulled from the
    # true first break to a later, geometrically smooth candidate. Keep each
    # breakpoint provisional and ask the next five usable traces which trend
    # their waveform candidates actually favour. Rejected breaks are repaired
    # one window at a time, so their error cannot propagate farther down-profile.
    maximum_confirmation_passes = max(1, len(order))
    for _ in range(maximum_confirmation_passes):
        y = np.asarray([selected[i]['index'] for i in order], dtype=float)
        prediction, breakpoints = _piecewise_linear_prediction(
            x, y, base_tolerance=base_tolerance,
            min_segment_length=min_segment)
        repaired = False
        previous_break = 0
        for breakpoint in breakpoints:
            trace_id = order[breakpoint]
            breakpoint_candidates = candidates_by_trace[trace_id]
            if (len(breakpoint_candidates) == 1 and
                    breakpoint_candidates[0].get('cost', 0.0) < 0):
                previous_break = breakpoint
                continue  # A user/validated anchor explicitly permits the break.

            history_start = max(previous_break, breakpoint-6)
            if breakpoint-history_start < 3:
                previous_break = breakpoint
                continue
            check_stop = min(len(order), breakpoint+6)
            check_positions = np.arange(breakpoint, check_stop, dtype=int)
            old_prediction = _robust_line_prediction(
                x[history_start:breakpoint], y[history_start:breakpoint],
                x[check_positions])

            old_support = 0
            new_support = 0
            # The triggering trace may simply be missing. Confirmation must
            # come from subsequent receivers, not from that trace alone.
            for local_id, position in enumerate(check_positions[1:], start=1):
                following_id = order[position]
                candidates = candidates_by_trace[following_id]
                old_candidate, old_score = _candidate_for_trend(
                    candidates, old_prediction[local_id], base_tolerance)
                new_candidate, new_score = _candidate_for_trend(
                    candidates, prediction[position], base_tolerance)
                if old_candidate is None and new_candidate is not None:
                    new_support += 1
                elif new_candidate is None and old_candidate is not None:
                    old_support += 1
                elif old_candidate is not None and new_candidate is not None:
                    if old_candidate['index'] == new_candidate['index']:
                        continue
                    if old_score+0.03 < new_score:
                        old_support += 1
                    elif new_score+0.03 < old_score:
                        new_support += 1

            if new_support >= 3 and new_support > old_support:
                previous_break = breakpoint
                continue

            changed = False
            for local_id, position in enumerate(check_positions):
                repair_id = order[position]
                candidates = candidates_by_trace[repair_id]
                if (len(candidates) == 1 and
                        candidates[0].get('cost', 0.0) < 0):
                    continue
                replacement, _ = _candidate_for_trend(
                    candidates, old_prediction[local_id], base_tolerance)
                forward_prediction = forward_predictions.get(repair_id)
                if (replacement is not None and
                        forward_prediction is not None and
                        abs(replacement['index']-forward_prediction) >
                        2.5*base_tolerance):
                    continue
                if (replacement is not None and
                        replacement['index'] != selected[repair_id]['index']):
                    selected[repair_id] = replacement
                    changed = True
            if changed:
                repaired = True
                break  # Refit globally before examining another breakpoint.
            previous_break = breakpoint
        if not repaired:
            break

    y = np.asarray([selected[i]['index'] for i in order], dtype=float)
    prediction, breakpoints = _piecewise_linear_prediction(
        x, y, base_tolerance=base_tolerance,
        min_segment_length=min_segment)
    segment_id = 0
    breakpoint_set = set(breakpoints)
    for position, trace_id in enumerate(order):
        if position in breakpoint_set:
            segment_id += 1
        candidate = dict(selected[trace_id])
        candidate_prediction = (forward_predictions[trace_id]
                                if trace_id in suspected_traces
                                else prediction[position])
        residual = abs(candidate['index']-candidate_prediction)
        candidate['spatial_prediction'] = float(candidate_prediction)
        candidate['spatial_residual'] = float(residual)
        candidate['spatial_agreement'] = float(np.exp(
            -residual/max(base_tolerance, np.finfo(float).eps)))
        candidate['segment_id'] = segment_id
        selected[trace_id] = candidate
    selected.update(independent)
    return selected

# Need to take a closer look at this: https://programmerall.com/article/10751929193/
# https://matplotlib.org/devdocs/gallery/widgets/polygon_selector_demo.html#polygon-selector (for the line selection)
# https://build-system.fman.io/ (for building into executable)


class MplCanvas(FigureCanvasQTAgg):
    """
    This class, MplCanvas, extends FigureCanvasQTAgg to create a custom canvas
    for matplotlib figures that can be integrated into a PyQt application.
    It sets up the figure, axes, and initial plot limits.
    """

    def __init__(self, parent=None, width=5, height=4, dpi=75):
        """Create a Qt-compatible Matplotlib canvas with one set of axes."""
        # Create a new figure with specified dimensions and resolution
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        # Add a single subplot to the figure
        self.axes = self.fig.add_subplot(111)
        # Initialize colorbar attribute (to be set later)
        self.cBar = None  # For the (eventual) colorbar
        # Initialize the FigureCanvasQTAgg with our figure
        super(MplCanvas, self).__init__(self.fig)
        # Set the parent widget for this canvas
        self.setParent(parent)
        # Define the default x and y limits for the plot
        self.homeXLimits = (0, 3)
        self.homeYLimits = (0, 24)


class CustomHomeToolbar(NavigationToolbar2QT):
    """
    This class, CustomHomeToolbar, extends NavigationToolbar2QT to create a custom toolbar
    with a modified 'home' functionality. It's designed to work with a custom canvas
    that has predefined 'home' limits for x and y axes.
    """

    def __init__(self, canvas, parent):
        """Attach the toolbar to a canvas that defines custom home limits."""
        # Initialize the parent class with the given canvas and parent
        super().__init__(canvas, parent)

    def home(self):
        """Restore the application-defined home limits and redraw the canvas."""
        # Override the default 'home' method
        # Set the x-axis limits to the predefined 'home' limits
        self.canvas.axes.set_xlim(*self.canvas.homeXLimits)
        # Set the y-axis limits to the predefined 'home' limits
        self.canvas.axes.set_ylim(*self.canvas.homeYLimits)
        # Redraw the canvas to reflect the changes
        self.canvas.draw_idle()

# class VerticalNavigationToolbar2QT(NavigationToolbar2QT):
#     def __init__(self, canvas, parent, coordinates=True):
#         super().__init__(canvas, parent, coordinates)


class paths:
    """
    Class for storing paths to files
    """

    def __init__(self) -> None:
        """
        Initialize the class instance variables
        """
        self.directory = []     # directory for the data files
        self.geometryFile = []  # geometry file path
        self.nbFiles = []       # number of SEG-2 files referenced in the geometry file
        self.seg2Files = []     # list of SEG-2 files


class geometry:
    """
    Class for storing geometry data
    """

    def __init__(self) -> None:
        """
        Initialize the class instance variables
        """
        self.sensors = []   # List of sensors
        self.sourcesId = []  # List of sources IDs
        self.receivers = []  # List of receivers


class model:
    """
    Class for storing inversion results
    """

    def __init__(self) -> None:
        """
        Initialize the class instance variables
        """
        self.nbLayers = 2
        self.thickLeft = np.ones((self.nbLayers-1,))*5
        self.dipAngles = np.zeros((self.nbLayers-1,))
        self.Vp = np.linspace(1000, 3000, num=self.nbLayers)

    def changeLayers(self, nbLayers):
        """
        Initialize the model for a given number of layers
        """
        self.nbLayers = nbLayers
        self.thickLeft = np.ones((self.nbLayers-1,))*5
        self.dipAngles = np.zeros((self.nbLayers-1,))
        self.Vp = np.linspace(1000, 3000, num=self.nbLayers)


class animationPicking():
    """Hold transient mouse and animation state for the picking interface."""

    def __init__(self) -> None:
        """Initialize the picker interaction state with neutral defaults."""
        self.timeOnClick = 0        # To know if the click is in fixed position of to pan/zoom
        # Storing the current position of the mouse
        self.mousePosition = [0, 0]
        self.currSelect = 0         # Storing the current trace number beiing analyzed
        self.changedSelect = True   # Variable to tell if the trace selected is different
        self.first = True           # Variable to tell if plotting for the first time
        # Maximum time (in sec) to consider a click on place
        self.maxClickLength = 0.5
        # Checking if currently picking of not (zooming or panning)
        self.notPicking = True
        # Checking that the fft is beiing displayed (true) or not (false)
        self.fftShowed = False
        self.fftAnim = None
        # Storing the initial position of the mouse when registering a mouse-click
        self.mousePositionInit = [0, 0]
        # Centre of the last time window selected by a click. It is kept
        # separately from mousePosition because the latter follows the cursor.
        self.autoPickCenter = None
        # Multi-pick selection on the currently displayed shot.
        self.selectedPicks = set()
        # Ctrl held when the current left-button drag started.
        self.linePicking = False


class inversionData():
    """Store inversion parameters, pyGIMLi objects, and the starting model."""

    def __init__(self) -> None:
        """Initialize default regularization, velocity, and mesh settings."""
        self.lam = 20.0
        self.zWeight = 0.5
        self.vTop = 500.0
        self.vBottom = 3000.0
        self.vMin = 10.0
        self.vMax = 5000.0
        self.startModel = None
        self.meshMaxCellSize = 5.0
        self.meshDepthMax = 50.0
        self.meshQuality = 33.33
        self.secNodes = 3
        self.blockyModel = False
        # Pygimli inversion features:
        self.mesh = None
        self.data = None
        self.manager = None

    def setStartModelGradient(self, data, mesh):
        """Create a vertical-gradient starting model on *mesh* for *data*."""
        self.startModel = pg.Vector(createGradientModel2D(
            data, mesh, self.vTop, self.vBottom))


class modellingAnimation():
    """Hold interaction state for editing the layered-model hodograph."""

    def __init__(self) -> None:
        """Initialize model-editing tolerances and selection state."""
        self.currPosition = [0, 0]
        self.maxClickLength = 0.5
        self.timeOnClick = 0
        self.offsetVelocity = 0.5
        self.tol = [1, 0.002]  # x in meteres, y in seconds
        self.currPointId = 1
        self.changingPts = False
        self.namesSources = []
        self.namesOrientations = []


class modellingData():
    """Store hodograph measurements and derived layered-model properties."""

    def __init__(self) -> None:
        """Initialize empty modelling collections and a two-layer default."""
        self.sensors = []
        self.measurements = []
        self.hodoPoints = []
        self.combinationSR = []
        self.model = []
        self.nbLayers = 2
        self.sourcesX = []
        self.orientations = []
        self.appVelocities = []
        self.interceptTime = []


class dataStorage():
    """Aggregate all mutable data and UI state owned by the main window."""

    def __init__(self) -> None:
        """Create empty acquisition, picking, modelling, and inversion state."""
        # Data variables:
        self.paths = paths()
        self.geometry = geometry()
        self.sisDataOriginal = []
        self.sisData = []
        self.picking = []
        self.pickingError = []
        self.manualPickingMask = []
        self.defectiveGeophones = []
        self.defectiveGeophoneReasons = []
        # Per-shot user override: 0 = automatic decision, 1 = defective,
        # -1 = explicitly restored as valid.
        self.manualDefectiveOverrides = []
        self.zeroOffsetTraces = []
        self.model = model()
        self.invData = inversionData()
        self.modellingAnimation = modellingAnimation()
        self.modellingData = modellingData()
        # Status variables:
        self.dataLoaded = False
        self.pickingDone = False
        self.inversionDone = False
        self.meshLoaded = False
        self.sisFileId = 0
        self.beginTime = []
        # Display-only multiplier for the wiggle amplitude in the gather.
        self.visualAmplitudeScale = 1.0
        # Graphical animation variables:
        self.animationPicking = animationPicking()


class picklingStatus():
    """Serializable subset of application state used by save/load actions."""

    def __init__(self) -> None:
        """Initialize the fields persisted in a picking-session file."""
        self.paths = paths()
        self.geometry = geometry()
        self.sisDataOriginal = []
        self.sisData = []
        self.beginTime = []
        self.picking = []
        self.pickingError = []
        self.manualPickingMask = []
        self.manualDefectiveOverrides = []
        self.sisFileId = 0


class PickT0(QDialog):
    '''
    This window intends to propose multiple options for the t0 picking.
        - Manual picking from the trace at the source (DONE)
        - Setting an offset value (DONE)
        - Using the picks that have been realise to interpolate the t0 at the source (DONE).
        - Disableing the functionalities for offset triggers.
        - Adapt zoom for pretrigger (zoom around the current t0 pick? - option ?).
    '''

    def __init__(self, parent):
        """Build the t0-correction dialog for the supplied main window."""
        super().__init__()
        self.setWindowTitle('Picking t0 helper')
        self.setWindowIcon(QtGui.QIcon(
            './images/SardineRebornLogo_100ppp.png'))
        self.resize(600, 300)
        self.mainWindow = parent  # In order to be able to plot elements and retreive values
        # From 0 to 999 negative values, 1000=0 from 1001 to 2000 positive
        self.nbValuesDisp = 2000
        self.initialUpdate = False
        # Create the file selector:
        # Comb Box for choosing the correct file to pick.
        self.comboBoxFilesPicking = QComboBox()
        self.sisFileId = self.mainWindow.dataUI.sisFileId
        for i, name in enumerate(self.mainWindow.dataUI.paths.seg2Files):
            self.comboBoxFilesPicking.addItem(name)
            if i == self.sisFileId:
                textItem = name
        self.comboBoxFilesPicking.setCurrentText(textItem)
        self.comboBoxFilesPicking.currentIndexChanged.connect(
            self.comboBoxChange)
        self.newT0 = np.zeros_like(
            self.mainWindow.dataUI.paths.seg2Files, dtype=float)

        # Creating the graph widget:
        self.traceGraph = MplCanvas(self, width=6, height=2)

        # Create the slider for graph:
        textZoom = QLabel('Zoom : ')
        self.sliderZoom = QSlider(QtCore.Qt.Horizontal, self)
        self.sliderZoom.setRange(1, 1000)
        self.sliderZoom.setValue(100)
        self.sliderZoom.valueChanged.connect(self.updateZoom)
        zoomLayout = QHBoxLayout()
        zoomLayout.addWidget(textZoom)
        zoomLayout.addWidget(self.sliderZoom)

        textPick = QLabel('T0 offset : ')
        self.slider = QSlider(QtCore.Qt.Horizontal, self)
        self.slider.setRange(0, self.nbValuesDisp)
        self.slider.setValue(int(self.nbValuesDisp/2))  # Initialize to "0.0"
        # self.slider.sliderReleased.connect(self.updateSlider)
        self.slider.valueChanged.connect(self.updateAxisSlider)
        textRange = QLabel('Range : ')
        self.rangeSlider = QDoubleSpinBox()
        self.rangeSlider.setDecimals(1)
        self.rangeSlider.setRange(0.1, 1.0)
        self.rangeSlider.setSingleStep(0.1)
        self.rangeSlider.setValue(0.1)
        self.rangeSlider.valueChanged.connect(self.updateRangeSlider)
        sliderLayout = QHBoxLayout()
        sliderLayout.addWidget(textPick)
        sliderLayout.addWidget(self.slider)
        sliderLayout.addWidget(textRange)
        sliderLayout.addWidget(self.rangeSlider)
        # Creating the widget with the current value (can be changed by the user):
        self.textLabel1 = QLabel('Value of time offset : ')
        self.currValue = QLineEdit(str(0), self)
        self.currValue.textEdited.connect(self.updateText)
        valueLayout = QHBoxLayout()
        valueLayout.addWidget(self.textLabel1)
        valueLayout.addWidget(self.currValue)

        # Creating the automated picking option:
        self.automatedPick = QPushButton('Automated', self)
        self.textLabel2 = QLabel('Number of Traces : ', self)
        self.nbTraces = QLineEdit(str(3), self)
        automatedLayout = QHBoxLayout()
        automatedLayout.addWidget(self.automatedPick)
        automatedLayout.addWidget(self.textLabel2)
        automatedLayout.addWidget(self.nbTraces)
        self.automatedPick.clicked.connect(self.automatedPicking)

        # Text line for update:
        self.textConfirm = QLabel(
            'To confirm selection, close the window.', self)
        self.textConfirm.setStyleSheet('background-color: cyan')
        confirmLayout = QHBoxLayout()
        confirmLayout.addWidget(self.textConfirm)
        confirmLayout.setAlignment(QtCore.Qt.AlignCenter)

        # Creating final layout:
        layout = QVBoxLayout()
        layout.addWidget(self.comboBoxFilesPicking)
        layout.addWidget(self.traceGraph)
        layout.addLayout(zoomLayout)
        layout.addLayout(sliderLayout)
        layout.addLayout(valueLayout)
        layout.addLayout(automatedLayout)
        layout.addLayout(confirmLayout)

        self.setLayout(layout)
        self.graphUpdate()

    def automatedPicking(self):
        '''
        This picking of the 0-time offset is done using the already picked traces in the dataset.
        It uses the closest n traces in any given directions to infer an hodochronique and 
        thus a 0-time offset.

        It is the only approach that might work for offset shots.
        '''
        sensors = self.mainWindow.dataUI.geometry.sensors
        sourcesId = self.mainWindow.dataUI.geometry.sourcesId
        receivers = self.mainWindow.dataUI.geometry.receivers
        currFile = self.sisFileId
        sId = int(sourcesId[currFile])
        distReceivers = calculateDistance(receivers, sensors[sId])
        k = int(self.nbTraces.text())
        pickedTimes = self.mainWindow.dataUI.picking[currFile, :]
        if min(distReceivers) == 0:
            idMin = np.argmin(distReceivers)
            distReceivers = np.delete(distReceivers, idMin)
            pickedTimes = np.delete(pickedTimes, idMin)
        idx = np.argpartition(distReceivers, k)[:k]
        # Retreive the picked times for the given files:
        pickedTime = pickedTimes[idx]
        if np.isnan(pickedTime).any():
            QMessageBox.warning(
                self, 'Warning !', 'Not all closests traces are already picked.\nUnable to infer the offset')
        else:
            _, timeAt0, _, _, _ = stats.linregress(
                distReceivers[idx], pickedTime)
            self.currValue.setText(str(timeAt0))
            self.newT0[self.sisFileId] = timeAt0
            if self.signal is not None:
                sliderPos = int((timeAt0-self.timeSEG2[0])/(
                    ((self.timeSEG2[-1]-self.timeSEG2[0])/1000*self.sliderZoom.value())-self.timeSEG2[0])*self.nbValuesDisp)
                self.slider.setValue(sliderPos)

    def graphUpdate(self):
        """Refresh the source-coincident trace and controls for the active shot."""
        # Changing the graph values:
        self.signal = None
        axTrace = self.traceGraph.axes
        axTrace.clear()
        sensors = self.mainWindow.dataUI.geometry.sensors
        sourcesId = self.mainWindow.dataUI.geometry.sourcesId
        receivers = self.mainWindow.dataUI.geometry.receivers
        currFile = self.sisFileId
        deltaT = float(self.mainWindow.dataUI.sisData[currFile][0].stats.delta)
        nbPoints = self.mainWindow.dataUI.sisData[currFile][0].stats.npts
        timeSEG2 = np.arange(
            self.mainWindow.dataUI.beginTime[currFile], self.mainWindow.dataUI.beginTime[currFile]+nbPoints*deltaT, deltaT)
        sId = int(sourcesId[currFile])
        found = False
        for i in range(len(self.mainWindow.dataUI.sisData[currFile])):
            rId = int(sensors.index(receivers[i]))
            if sId == rId:
                found = True
                self.deltaT = deltaT
                self.signal = self.mainWindow.dataUI.sisData[currFile][i]
                self.timeSEG2 = timeSEG2
                break
        if found:
            # The source is located at the same position as a single trace
            self.sliderZoom.setEnabled(True)
            self.slider.setEnabled(True)
            axTrace.plot(timeSEG2, self.signal)
            if self.newT0[self.sisFileId] == 0:
                currPicking = self.mainWindow.dataUI.picking[currFile, i]
                if not (np.isnan(currPicking)):
                    axTrace.axvline(currPicking, color='g')
                    rangeValue = self.rangeSlider.value()
                    sliderPos = int(currPicking/rangeValue *
                                    self.nbValuesDisp/2 + self.nbValuesDisp/2)
                    if sliderPos < 0:
                        sliderPos = 0
                        currPicking = -rangeValue
                    elif sliderPos > self.nbValuesDisp:
                        sliderPos = self.nbValuesDisp
                        currPicking = rangeValue
                    self.slider.setValue(sliderPos)
                    self.currValue.setText(str(currPicking))
                else:
                    axTrace.axvline(0, color='g')
                    currPicking = 0
                    sliderPos = int(self.nbValuesDisp/2)
                    self.initialUpdate = True
                    self.slider.setValue(sliderPos)
                    self.initialUpdate = False
            else:
                currPicking = self.newT0[self.sisFileId]
                axTrace.axvline(currPicking, color='g')
                rangeValue = self.rangeSlider.value()
                sliderPos = int(currPicking/rangeValue *
                                self.nbValuesDisp/2 + self.nbValuesDisp/2)
                if sliderPos < 0:
                    sliderPos = 0
                    currPicking = -rangeValue
                elif sliderPos > self.nbValuesDisp:
                    sliderPos = self.nbValuesDisp
                    currPicking = rangeValue
                self.slider.setValue(sliderPos)
                self.currValue.setText(str(currPicking))
            # axTrace.set_xlim(self.timeSEG2[0], (self.timeSEG2[-1]-self.timeSEG2[0])/1000*self.sliderZoom.value())
            axTrace.set_xlim(currPicking - self.sliderZoom.value() /
                             1000, currPicking + self.sliderZoom.value()/1000)
        else:
            # The source is not located at the position of a receiver. t0 cannot be picked from the graph.
            self.sliderZoom.setEnabled(False)
            self.slider.setEnabled(False)
        self.traceGraph.draw()

    def updateZoom(self):
        """Center the trace view on the current t0 value at the chosen scale."""
        axTrace = self.traceGraph.axes
        currPicking = float(self.currValue.text())
        # self.timeSEG2[0], (self.timeSEG2[-1]-self.timeSEG2[0])/1000*self.sliderZoom.value())
        axTrace.set_xlim(currPicking - self.sliderZoom.value() /
                         1000, currPicking + self.sliderZoom.value()/1000)
        self.traceGraph.draw()
        # # Recompute value slider position:
        # currPick = float(self.currValue.text())
        # sliderPos = int((currPick-self.timeSEG2[0])/(((self.timeSEG2[-1]-self.timeSEG2[0])/1000*self.sliderZoom.value())-self.timeSEG2[0])*self.nbValuesDisp)
        # self.initialUpdate = True
        # self.slider.setValue(sliderPos)
        # self.initialUpdate = False

    def comboBoxChange(self, newId):
        """Select another shot by index and refresh its t0 preview."""
        self.sisFileId = newId
        self.graphUpdate()

    # def updateSlider(self):
    #     pass
    #     # currSlider = self.slider.value()
    #     # currPick = self.timeSEG2[0] + currSlider*(((self.timeSEG2[-1]-self.timeSEG2[0])/1000*self.sliderZoom.value())-self.timeSEG2[0])/self.nbValuesDisp
    #     # self.currValue.setText(str(currPick))
    #     # self.newT0[self.sisFileId] = float(currPick)

    def updateAxisSlider(self):
        """Convert the t0 slider position to seconds and redraw the marker."""
        if self.signal is not None:
            axTrace = self.traceGraph.axes
            # The source is located at the same position as a single trace
            axTrace.clear()
            currSlider = self.slider.value()
            rangeValue = self.rangeSlider.value()
            currPicking = ((currSlider - self.nbValuesDisp/2) /
                           (self.nbValuesDisp/2))*rangeValue
            axTrace.plot(self.timeSEG2, self.signal)
            axTrace.axvline(currPicking, color='g')
            # self.timeSEG2[0], (self.timeSEG2[-1]-self.timeSEG2[0])/1000*self.sliderZoom.value())
            axTrace.set_xlim(currPicking - self.sliderZoom.value() /
                             1000, currPicking + self.sliderZoom.value()/1000)
            # axTrace.set_xlim(self.timeSEG2[0], (self.timeSEG2[-1]-self.timeSEG2[0])/1000*self.sliderZoom.value())
            self.traceGraph.draw()
            # Change values in the text box and back-end:
            if not (self.initialUpdate):
                self.currValue.setText(str(currPicking))
                self.newT0[self.sisFileId] = float(currPicking)

    def updateText(self):
        """Validate a typed t0 value and synchronize the position slider."""
        try:
            currPicking = float(self.currValue.text())
            rangeValue = self.rangeSlider.value()
            sliderPos = int(currPicking/rangeValue *
                            self.nbValuesDisp/2 + self.nbValuesDisp/2)
            if sliderPos < 0:
                sliderPos = 0
                currPicking = -rangeValue
            elif sliderPos > self.nbValuesDisp:
                sliderPos = self.nbValuesDisp
                currPicking = rangeValue
            self.slider.setValue(sliderPos)
            self.newT0[self.sisFileId] = float(currPicking)
        except:
            pass

    def updateRangeSlider(self):
        """Rescale the t0 position slider while preserving the current value."""
        rangeValue = self.rangeSlider.value()  # New range value
        # We need to update the slider min/max values accordingly.
        currPicking = float(self.currValue.text())
        sliderPos = int(currPicking/rangeValue *
                        self.nbValuesDisp/2 + self.nbValuesDisp/2)
        if sliderPos < 0:
            sliderPos = 0
        elif sliderPos > self.nbValuesDisp:
            sliderPos = self.nbValuesDisp
        self.initialUpdate = True
        self.slider.setValue(sliderPos)
        self.initialUpdate = False

    def getNewT0(self):
        """Return the per-shot t0 corrections after a consistency check."""
        # Check that newT0 is up to date (necessary ?)
        if self.signal is not None:
            currSlider = self.slider.value()
            currPickSlider = self.timeSEG2[0] + currSlider*(
                ((self.timeSEG2[-1]-self.timeSEG2[0])/1000*self.sliderZoom.value())-self.timeSEG2[0])/self.nbValuesDisp
            sliderPrecision = (self.timeSEG2[0] + (currSlider+1)*(((self.timeSEG2[-1]-self.timeSEG2[0]) /
                               1000*self.sliderZoom.value())-self.timeSEG2[0])/self.nbValuesDisp) - currPickSlider
            currPickText = float(self.currValue.text())
            if abs(currPickSlider - currPickText) < sliderPrecision:
                if abs(self.newT0[self.sisFileId] - (currPickSlider+currPickText)/2) < sliderPrecision:
                    pass
                else:
                    QMessageBox.warning(
                        self, 'Warning!', 'Possible error while executing.\nCheck the sismograms.')
            else:
                QMessageBox.warning(
                    self, 'Warning!', 'Possible error while executing.\nCheck the sismograms.')
        return self.newT0


class Window(QMainWindow):
    """Main Sardine Reborn window and controller for the processing workflow."""

    def __init__(self) -> None:
        """Build the menus, tabs, shared state, plots, and event connections."""
        super().__init__()

        # Initializing the data structure
        self.dataUI = dataStorage()

        # Initialize the UI
        self.setWindowTitle('Sardine Reborn')
        self.setWindowIcon(QtGui.QIcon(
            './images/SardineRebornLogo_100ppp.png'))
        self.resize(1100, 600)

        # Adding tabs to the layout
        self.tabs = QTabWidget(self)
        self.tabs.addTab(self._pickTracesTabUI(), 'Traces picking')
        self.tabs.addTab(self._inversionTabUI(), 'Inversion')
        self.tabs.addTab(self._displayInversionResultsTabUI(),
                         'Display inversion results')
        self.tabs.addTab(self._modelTabUI(), 'Model')

        self.setCentralWidget(self.tabs)

        # Defining menu actions:
        # File openning / saving picking
        openFile = QAction("&Open Geometry File", self)
        openFile.setShortcut('Ctrl+O')  # Setting ctrl+o as the shortcut
        openFile.setStatusTip('Open the *.geometry file')
        openFile.triggered.connect(self._openGeometry)

        savePicking = QAction("Save Current &Picking", self)
        savePicking.setShortcut('Ctrl+P')
        savePicking.setStatusTip('Save the picking state in a *.sgt file')
        savePicking.triggered.connect(self._savePicking)

        loadPicking = QAction("&Load Picking File", self)
        loadPicking.setShortcut('Ctrl+L')
        loadPicking.setStatusTip('Load an existing *.sgt file for inversion')
        loadPicking.triggered.connect(self._loadPicking)

        savePickingStatus = QAction("Save the picking state", self)
        savePickingStatus.setStatusTip('Save the status into a pickle')
        savePickingStatus.triggered.connect(self.saveStatePicking)

        loadPickingStatus = QAction("Load a picking state", self)
        loadPickingStatus.setStatusTip('Load a pickle status')
        loadPickingStatus.triggered.connect(self.loadStatePicking)

        recalculateT0 = QAction("Recalculate t=0 from geometry", self)
        recalculateT0.setStatusTip(
            'Set t(distance=0)=0 from a co-located receiver or a robust extrapolation.')
        recalculateT0.triggered.connect(self.recalculateT0FromGeometry)

        # Filtering of the signal:
        fftGraph = QAction("Show FFT of dataset", self)
        fftGraph.setStatusTip('Show the FFT transform of the different traces')
        fftGraph.triggered.connect(self.fftGraph)

        dcFilter = QAction("Apply a DC filter", self)
        dcFilter.setStatusTip('Apply a DC filter to the signal')
        dcFilter.triggered.connect(self.dcFilter)

        trimFilter = QAction("Trim the datasets", self)
        trimFilter.setStatusTip('Apply trimming to the signal')
        trimFilter.triggered.connect(self.trimFilter)

        highpassFilter = QAction("Apply a high-pass filter", self)
        highpassFilter.setStatusTip('Apply a high pass filter to the signal')
        highpassFilter.triggered.connect(self.highpassFilter)

        lowpassFilter = QAction("Apply a low-pass filter", self)
        lowpassFilter.setStatusTip('Apply a low pass filter to the signal')
        lowpassFilter.triggered.connect(self.lowpassFilter)

        resetFilters = QAction("Reset the filters", self)
        resetFilters.setStatusTip('Reset all previously applied filters')
        resetFilters.triggered.connect(self.resetFilters)

        autoPicking = QAction("Automated picking", self)
        autoPicking.setStatusTip(
            'Hybrid STA/LTA + AIC first-break picking with multi-trace consistency.')
        autoPicking.triggered.connect(self.autoPicking)

        # Inversion through pigimli:
        loadInvMesh = QAction('Load Inversion Mesh', self)
        loadInvMesh.setStatusTip(
            'Load an inversion mesh that was already created (*.poly)')
        loadInvMesh.triggered.connect(self._loadInvMesh)

        loadInitialModel = QAction('Load Initial Model', self)
        loadInitialModel.setStatusTip(
            'Load an existing model as the starting model (*.vector) for the inversion')
        loadInitialModel.triggered.connect(self._loadInitModel)

        saveInvMesh = QAction('Save Inversion Mesh', self)
        saveInvMesh.setStatusTip('Save the inversion mesh (*.poly)')
        saveInvMesh.triggered.connect(self._saveInvMesh)

        saveInvAsVTK = QAction('Save Inversion as VTK', self)
        saveInvAsVTK.setStatusTip(
            'Save the inversion results into a VTK file (for Paraview)')
        saveInvAsVTK.triggered.connect(self._saveInvAsVTK)

        saveInvResponse = QAction('Save the Inverse Response', self)
        saveInvResponse.setStatusTip(
            'Save the model response for the last iteration (*.vector)')
        saveInvResponse.triggered.connect(self._saveInvResponse)

        saveInvResult = QAction('Save the Inverse Results', self)
        saveInvResult.setStatusTip(
            'Save the model for the last iteration (*.vector)')
        saveInvResult.triggered.connect(self._saveInvResult)

        # saveModel = QAction("Save Current &Model",self)
        # saveModel.setShortcut('Ctrl+M')
        # saveModel.setStatusTip('Save the current model in a *.txt file')
        # saveModel.triggered.connect(self._saveModel)

        # Adding the menu bar atop
        menuBarInternal = self.menuBar()
        menuBarInternal.setNativeMenuBar(False)
        fileMenu = menuBarInternal.addMenu("Picking")
        fileMenu.addAction(openFile)
        fileMenu.addAction(savePicking)
        fileMenu.addAction(loadPicking)
        fileMenu.addAction(savePickingStatus)
        fileMenu.addAction(loadPickingStatus)
        fileMenu.addAction(recalculateT0)
        fileMenu = menuBarInternal.addMenu("Filters")
        fileMenu.addAction(fftGraph)
        fileMenu.addAction(dcFilter)
        fileMenu.addAction(trimFilter)
        fileMenu.addAction(highpassFilter)
        fileMenu.addAction(lowpassFilter)
        fileMenu.addAction(resetFilters)
        fileMenu.addAction(autoPicking)
        # fileMenu.addAction(saveModel)
        fileMenu = menuBarInternal.addMenu("Inversion")
        fileMenu.addAction(loadInvMesh)
        fileMenu.addAction(loadInitialModel)
        fileMenu.addAction(saveInvMesh)
        fileMenu.addAction(saveInvResult)
        fileMenu.addAction(saveInvResponse)
        fileMenu.addAction(saveInvAsVTK)

        # Defining the status bar:
        self.statusBar = QStatusBar(self)
        permanentMessage = QLabel(self.statusBar)
        permanentMessage.setText('Sardine Reborn - Hadrien Michel (2023)')
        self.statusBar.addPermanentWidget(permanentMessage)
        self.setStatusBar(self.statusBar)
        # Adding easy access to fftGraph action for changes in name and tip
        self.dataUI.animationPicking.fftAnim = fftGraph
        # --- Fine tuning picking with keyboard arrows ---
        self.sc_pick_left = QShortcut(QKeySequence(QtCore.Qt.Key_Left), self)
        self.sc_pick_right = QShortcut(QKeySequence(QtCore.Qt.Key_Right), self)

        # Important: capture even if focus is in children (tabs, widgets, canvas…)
        self.sc_pick_left.setContext(QtCore.Qt.ApplicationShortcut)
        self.sc_pick_right.setContext(QtCore.Qt.ApplicationShortcut)

        self.sc_pick_left.activated.connect(lambda: self.nudge_pick(-1))
        self.sc_pick_right.activated.connect(lambda: self.nudge_pick(+1))

        self.sc_pick_left_fast = QShortcut(QKeySequence("Shift+Left"), self)
        self.sc_pick_right_fast = QShortcut(QKeySequence("Shift+Right"), self)
        self.sc_pick_left_fast.setContext(QtCore.Qt.ApplicationShortcut)
        self.sc_pick_right_fast.setContext(QtCore.Qt.ApplicationShortcut)
        self.sc_pick_left_fast.activated.connect(lambda: self.nudge_pick(-1, n_samples=10))
        self.sc_pick_right_fast.activated.connect(lambda: self.nudge_pick(+1, n_samples=10))
        QApplication.instance().installEventFilter(self)

    def eventFilter(self, obj, event):
        """
        Intercepts keyboard events *before* they reach Qt widgets.

        Purpose:
        - LEFT / RIGHT fine-tune seismic time picking
        - UP / DOWN navigate through traces
        - A automatically picks inside the currently displayed zoom window
        - DELETE removes the current pick
        - Works even if a QSpinBox / QLineEdit ("Trace Number") has the focus
        - Active only in the 'Traces picking' tab
        """
        if event.type() == QtCore.QEvent.KeyPress:
            key = event.key()
            if self.tabs.currentIndex() == 0:
                if key in (QtCore.Qt.Key_Left, QtCore.Qt.Key_Right):
                    # FINE STEP: arrow key only
                    if event.modifiers() == QtCore.Qt.NoModifier:
                        direction = -1 if key == QtCore.Qt.Key_Left else +1
                        self.nudge_pick(direction, n_samples=1)
                        event.accept()
                        return True

                    # COARSE STEP: Shift + arrow
                    if event.modifiers() == QtCore.Qt.ShiftModifier:
                        direction = -1 if key == QtCore.Qt.Key_Left else +1
                        self.nudge_pick(direction, n_samples=10)
                        event.accept()
                        return True

                if key in (QtCore.Qt.Key_Up, QtCore.Qt.Key_Down) and event.modifiers() == QtCore.Qt.NoModifier:
                    # Trace 0 is drawn at the bottom of the gather: Up moves
                    # to the visually higher trace, Down to the lower trace.
                    direction = +1 if key == QtCore.Qt.Key_Up else -1
                    self.change_current_trace(direction)
                    event.accept()
                    return True

                if key == QtCore.Qt.Key_A and event.modifiers() == QtCore.Qt.NoModifier:
                    if not event.isAutoRepeat():
                        self.auto_pick_current_window()
                    event.accept()
                    return True

                if key in (QtCore.Qt.Key_Delete, QtCore.Qt.Key_Backspace):
                    self.delete_current_pick()
                    event.accept()
                    return True

        return super().eventFilter(obj, event)


    def nudge_pick(self, direction: int, n_samples: int = 1):
        """
        direction: -1 (gauche) ou +1 (droite)
        n_samples: nombre d'échantillons à déplacer (1 = ultra fin)
        """
        # Sécurité: pas de données
        if len(self.dataUI.sisData) == 0:
            return

        file_id = self.dataUI.sisFileId
        tr_id = self.dataUI.animationPicking.currSelect

        # Sécurité: index out of range
        if file_id is None or tr_id is None:
            return

        # Pick actuel
        curr = self.dataUI.picking[file_id, tr_id]
        if np.isnan(curr):
            return  # rien à affiner si pas encore pické

        # Pas de temps = 1 sample (super fin)
        dt = float(self.dataUI.sisData[file_id][0].stats.delta)  # seconds
        step = direction * n_samples * dt

        new_pick = curr + step

        # Empêcher d'aller avant beginTime (chez toi beginTime sert de "t0"/offset)
        min_t = float(self.dataUI.beginTime[file_id])
        if new_pick < min_t:
            new_pick = min_t

        self.ensure_manual_picking_mask()
        self.dataUI.picking[file_id, tr_id] = new_pick
        self.dataUI.manualPickingMask[file_id, tr_id] = True


        # (Optionnel) tu peux aussi ajuster l’erreur si tu veux, sinon ne touche pas
        # self.dataUI.pickingError[file_id, tr_id] = self.dataUI.pickingError[file_id, tr_id]

        # Forcer le refresh des graphs via ton système d’animation
        self.dataUI.animationPicking.changedSelect = True

    def ensure_manual_picking_mask(self):
        """Create a provenance mask, treating legacy existing picks as manual."""
        mask = np.asarray(self.dataUI.manualPickingMask)
        if mask.shape != np.asarray(self.dataUI.picking).shape:
            self.dataUI.manualPickingMask = np.isfinite(
                self.dataUI.picking).astype(bool)

    def ensure_manual_defective_overrides(self):
        """Create the per-shot defective-trace overrides for legacy states."""
        shape = np.asarray(self.dataUI.picking).shape
        overrides = np.asarray(self.dataUI.manualDefectiveOverrides)
        if overrides.shape != shape:
            self.dataUI.manualDefectiveOverrides = np.zeros(shape, dtype=np.int8)

    def effective_defective_mask(self, file_id):
        """Combine automatic global detection with per-shot user decisions."""
        self.ensure_manual_defective_overrides()
        nb_traces = len(self.dataUI.sisData[file_id])
        automatic = np.asarray(self.dataUI.defectiveGeophones, dtype=bool)
        if automatic.shape != (nb_traces,):
            automatic = np.zeros(nb_traces, dtype=bool)
        overrides = self.dataUI.manualDefectiveOverrides[file_id]
        return np.where(overrides == 1, True,
                        np.where(overrides == -1, False, automatic))

    def defective_trace_reason(self, file_id, trace_id):
        """Return the display reason for an effective defective trace."""
        self.ensure_manual_defective_overrides()
        override = self.dataUI.manualDefectiveOverrides[file_id, trace_id]
        if override == 1:
            return 'manually marked defective'
        if len(self.dataUI.defectiveGeophoneReasons) > trace_id:
            return self.dataUI.defectiveGeophoneReasons[trace_id]
        return ''

    def manage_current_defective_trace(self):
        """Mark or restore the selected trace for this shot or every shot."""
        if not self.dataUI.dataLoaded or len(self.dataUI.sisData) == 0:
            return
        self.ensure_manual_defective_overrides()
        trace_id = self.dataUI.animationPicking.currSelect
        file_id = self.dataUI.sisFileId
        choices = [
            'Mark defective for this shot',
            'Mark defective for all shots',
            'Restore as valid for this shot',
            'Restore as valid for all shots',
        ]
        choice, ok = QInputDialog.getItem(
            self, 'Trace validity',
            f'Geophone {trace_id+1}: choose its validity status.',
            choices, 0, False)
        if not ok:
            return
        is_all_shots = choice in (choices[1], choices[3])
        value = 1 if choice in (choices[0], choices[1]) else -1
        target_files = (np.arange(len(self.dataUI.sisData), dtype=int)
                        if is_all_shots else np.asarray([file_id], dtype=int))
        self.dataUI.manualDefectiveOverrides[target_files, trace_id] = value
        if value == 1:
            self.dataUI.picking[target_files, trace_id] = np.nan
            self.dataUI.pickingError[target_files, trace_id] = np.nan
            self.ensure_manual_picking_mask()
            self.dataUI.manualPickingMask[target_files, trace_id] = False
            action = 'marked defective'
        else:
            action = 'restored as valid'
        self.dataUI.animationPicking.selectedPicks.discard(trace_id)
        self.dataUI.animationPicking.changedSelect = True
        scope = 'all shots' if is_all_shots else 'this shot'
        self.statusBar.showMessage(
            f'Geophone {trace_id+1} {action} for {scope}.', 6000)

    def delete_current_pick(self):
        """Remove selected picks, or the current pick when no group is selected."""
        if len(self.dataUI.sisData) == 0:
            return
        file_id = self.dataUI.sisFileId
        self.ensure_manual_picking_mask()
        selected = self.dataUI.animationPicking.selectedPicks
        trace_ids = (sorted(selected) if selected else
                     [self.dataUI.animationPicking.currSelect])
        for trace_id in trace_ids:
            self.dataUI.picking[file_id, trace_id] = np.nan
            self.dataUI.pickingError[file_id, trace_id] = np.nan
            self.dataUI.manualPickingMask[file_id, trace_id] = False
        selected.clear()
        self.dataUI.animationPicking.autoPickCenter = None
        self.dataUI.animationPicking.changedSelect = True
        self.statusBar.showMessage(
            f'{len(trace_ids)} pick(s) removed.', 3000)

    def select_picks_in_rectangle(self, x_values, y_values):
        """Select existing picks whose centres lie in the dragged rectangle."""
        x_min, x_max = sorted(x_values)
        y_min, y_max = sorted(y_values)
        picks = self.dataUI.picking[self.dataUI.sisFileId]
        selected = {
            trace_id for trace_id, pick in enumerate(picks)
            if np.isfinite(pick) and x_min <= pick <= x_max and
            y_min-0.5 <= trace_id <= y_max+0.5}
        self.dataUI.animationPicking.selectedPicks = selected
        self.dataUI.animationPicking.changedSelect = True
        self.statusBar.showMessage(
            f'{len(selected)} pick(s) selected. Press Delete to remove them.',
            5000)

    def change_current_trace(self, direction: int):
        """Move to the previous/next trace while remaining in picking mode."""
        if len(self.dataUI.sisData) == 0 or not self.spinBoxCurrSelect.isEnabled():
            return
        new_id = np.clip(
            self.dataUI.animationPicking.currSelect + direction,
            self.spinBoxCurrSelect.minimum(),
            self.spinBoxCurrSelect.maximum())
        self.spinBoxCurrSelect.setValue(int(new_id))

    def current_signal_time_bounds(self):
        """Return the complete time range of the currently selected file."""
        file_id = self.dataUI.sisFileId
        trace = self.dataUI.sisData[file_id][0]
        begin = float(self.dataUI.beginTime[file_id])
        dt = float(trace.stats.delta)
        return begin, begin + (trace.stats.npts - 1) * dt, dt

    def configure_time_window(self, reset=False):
        """Adapt the global-display controls to the current signal bounds."""
        if len(self.dataUI.sisData) == 0:
            return
        signal_start, signal_end, dt = self.current_signal_time_bounds()
        old_start = self.timeWindowStart.value()
        old_end = self.timeWindowEnd.value()

        for spinbox in (self.timeWindowStart, self.timeWindowEnd):
            spinbox.blockSignals(True)
            spinbox.setRange(signal_start, signal_end)
            spinbox.setSingleStep(max(dt, 10.0 * dt))

        if reset or old_end <= old_start:
            window_start, window_end = signal_start, signal_end
        else:
            window_start = np.clip(old_start, signal_start, signal_end)
            window_end = np.clip(old_end, signal_start, signal_end)
            if window_end - window_start < dt:
                window_start, window_end = signal_start, signal_end

        self.timeWindowStart.setValue(float(window_start))
        self.timeWindowEnd.setValue(float(window_end))
        self.timeWindowStart.blockSignals(False)
        self.timeWindowEnd.blockSignals(False)
        self.timeWindowStart.setEnabled(True)
        self.timeWindowEnd.setEnabled(True)
        self.buttonApplyTimeWindow.setEnabled(True)
        self.buttonResetTimeWindow.setEnabled(True)
        self.apply_time_window()

    def apply_time_window(self):
        """Apply the requested time limits to the global trace display."""
        if len(self.dataUI.sisData) == 0:
            return
        window_start = self.timeWindowStart.value()
        window_end = self.timeWindowEnd.value()
        _, _, dt = self.current_signal_time_bounds()
        if window_end - window_start < dt:
            self.statusBar.showMessage(
                'The end of the time window must be greater than its start.', 4000)
            return
        self.mainGraph.homeXLimits = (window_start, window_end)
        self.mainGraph.axes.set_xlim(window_start, window_end)
        self.dataUI.animationPicking.changedSelect = True
        self.statusBar.showMessage(
            f'Displayed time window: {window_start:.6f} to {window_end:.6f} s.', 3000)

    def reset_time_window(self):
        """Restore the complete time range of the current signal."""
        self.configure_time_window(reset=True)

    def auto_pick_current_window(self):
        """Pick with AIC inside the zoom window currently being displayed."""
        if not AIC_SIMPLE_IMPORT:
            QMessageBox.warning(
                self, 'Warning!',
                'The automatic AIC picker is unavailable in this ObsPy installation.')
            return
        if len(self.dataUI.sisData) == 0:
            return

        file_id = self.dataUI.sisFileId
        trace_id = self.dataUI.animationPicking.currSelect
        # The zoom canvas is centred on mousePosition. Therefore A must use
        # that same position, whether or not the user clicked beforehand.
        center = self.dataUI.animationPicking.mousePosition[0]
        if center is None or not np.isfinite(center):
            self.statusBar.showMessage(
                'Move the cursor over the signal before pressing A.', 4000)
            return

        trace = self.dataUI.sisData[file_id][trace_id]
        dt = float(trace.stats.delta)
        begin = float(self.dataUI.beginTime[file_id])
        data = np.asarray(trace.data, dtype=float)

        center_idx = int(round((center - begin) / dt))
        start_idx = max(0, center_idx - 100)
        stop_idx = min(data.size, center_idx + 101)
        window = data[start_idx:stop_idx]
        if window.size < 8 or not np.isfinite(window).all():
            self.statusBar.showMessage(
                'The selected window is too short or contains invalid values.', 4000)
            return

        # Ignore edge minima, which are commonly AIC boundary artefacts.
        aic = np.asarray(aic_simple(window), dtype=float)
        margin = max(2, int(round(0.05 * aic.size)))
        valid = aic[margin:aic.size-margin]
        if valid.size == 0 or not np.isfinite(valid).any():
            self.statusBar.showMessage('No AIC minimum found in this window.', 4000)
            return
        local_idx = margin + int(np.nanargmin(valid))
        picked_time = begin + (start_idx + local_idx) * dt

        self.ensure_manual_picking_mask()
        self.dataUI.picking[file_id, trace_id] = picked_time
        self.dataUI.pickingError[file_id, trace_id] = DEFAULT_ERROR
        self.dataUI.manualPickingMask[file_id, trace_id] = True
        self.dataUI.animationPicking.autoPickCenter = picked_time
        self.dataUI.animationPicking.mousePosition[0] = picked_time
        self.dataUI.animationPicking.changedSelect = True
        self.statusBar.showMessage(
            f'Trace {trace_id}: automatic AIC pick at {picked_time:.6f} s.', 4000)

    # Oppening and closing message boxes:
    def showEvent(self, event):
        """Show the application welcome message when the window opens."""
        msgBox = QMessageBox(self)
        msgBox.setIconPixmap(QtGui.QPixmap('./images/SardineRebornLogo_100ppp.png').scaled(
            200, 100, aspectRatioMode=QtCore.Qt.KeepAspectRatio))
        msgBox.setText('Welcome to Sardine Reborn!')
        msgBox.setInformativeText('by Hadrien Michel (2022)')
        msgBox.setWindowTitle('Welcome!')
        msgBox.show()
        event.accept()

    def closeEvent(self, event):
        """Ask for confirmation and close Matplotlib figures before exiting."""
        reply = QMessageBox.question(self, 'Closing ...', 'Are you sure you want to quit?',
                                     QMessageBox.Ok | QMessageBox.Cancel, QMessageBox.Ok)
        if reply == QMessageBox.Ok:
            event.accept()
            # To get rid of the figure openning after the closing.
            pyplot.close('all')
        else:
            event.ignore()

    # Animations definitions:
    # PyQt5 animations

    # Matplotlib animations:
    def changeMouse(self, event):
        """Record the latest cursor coordinates over the picking axes."""
        if event.inaxes is not None:
            self.dataUI.animationPicking.mousePosition = [
                event.xdata, event.ydata]
        return 0

    def change_visual_amplitude(self, event):
        """Change only the displayed gather amplitude with the mouse wheel."""
        if (not self.dataUI.dataLoaded or
                self.dataUI.animationPicking.fftShowed or
                event.inaxes is not self.mainGraph.axes):
            return 0
        if event.button == 'up':
            factor = 1.15
        elif event.button == 'down':
            factor = 1.0/1.15
        else:
            return 0
        self.dataUI.visualAmplitudeScale = float(np.clip(
            self.dataUI.visualAmplitudeScale*factor, 0.10, 10.0))
        self.visualAmplitudeLabel.setText(
            f'Visual amplitude: ×{self.dataUI.visualAmplitudeScale:.2f} '
            '(mouse wheel over the gather)')
        self.dataUI.animationPicking.changedSelect = True
        return 0

    def onPress(self, event):  # For both windows
        """Start a click or drag interaction in the seismic picking view."""
        if event.button == MouseButton.LEFT or event.button == MouseButton.RIGHT:
            self.dataUI.animationPicking.timeOnClick = time.time()
            # Checking if the zoom or pan tools are checked to enable line-picking
            toolbarState = str(self.mainGraph.fig.canvas.toolbar.mode)
            # toolbarOptions = type(toolbarState)
            if toolbarState != '':  # toolbarOptions.ZOOM or toolbarState == toolbarOptions.PAN:
                self.dataUI.animationPicking.notPicking = True
            else:
                self.dataUI.animationPicking.notPicking = False
                self.dataUI.animationPicking.mousePositionInit = self.dataUI.animationPicking.mousePosition
                modifiers = QApplication.keyboardModifiers()
                self.dataUI.animationPicking.linePicking = bool(
                    event.button == MouseButton.LEFT and
                    modifiers & QtCore.Qt.ControlModifier)
        return 0

    def onRelease(self, event):
        """Commit a point, line, selection rectangle, or pick-error gesture."""
        if not (self.dataUI.animationPicking.notPicking):
            # If left click and not dragging accross the pannel
            if event.button == MouseButton.LEFT and ((time.time() - self.dataUI.animationPicking.timeOnClick) < self.dataUI.animationPicking.maxClickLength):
                self.dataUI.animationPicking.selectedPicks.clear()
                if self.buttonTabPickingSetT0.isChecked():
                    self.dataUI.beginTime[self.dataUI.sisFileId] -= self.dataUI.animationPicking.mousePosition[0]
                    # Change all picked signals for this stream:
                    self.dataUI.picking[self.dataUI.sisFileId,
                                        :] -= self.dataUI.animationPicking.mousePosition[0]
                    self.buttonTabPickingSetT0.setChecked(False)
                else:
                    # Reject only clicks before the recorded signal. Negative
                    # arrival times are valid before a t0 correction.
                    if self.dataUI.animationPicking.mousePosition[0] < self.dataUI.beginTime[self.dataUI.sisFileId]:
                        self.dataUI.picking[self.dataUI.sisFileId,
                                            self.dataUI.animationPicking.currSelect] = np.nan
                        self.dataUI.pickingError[self.dataUI.sisFileId,
                                                 self.dataUI.animationPicking.currSelect] = np.nan
                        self.ensure_manual_picking_mask()
                        self.dataUI.manualPickingMask[
                            self.dataUI.sisFileId,
                            self.dataUI.animationPicking.currSelect] = False
                        self.dataUI.animationPicking.autoPickCenter = None
                        QMessageBox.warning(
                            self, 'Warning!', 'The selected time is outside the recorded signal.')
                    else:
                        self.ensure_manual_picking_mask()
                        self.dataUI.picking[self.dataUI.sisFileId,
                                            self.dataUI.animationPicking.currSelect] = self.dataUI.animationPicking.mousePosition[0]
                        self.dataUI.animationPicking.autoPickCenter = self.dataUI.animationPicking.mousePosition[0]
                        # max(self.dataUI.animationPicking.mousePosition[0]*DEFAULT_ERROR, 0.000001) # Default error is 3%
                        self.dataUI.pickingError[self.dataUI.sisFileId,
                                                 self.dataUI.animationPicking.currSelect] = DEFAULT_ERROR
                        self.dataUI.manualPickingMask[
                            self.dataUI.sisFileId,
                            self.dataUI.animationPicking.currSelect] = True
                self.dataUI.animationPicking.changedSelect = True
            elif event.button == MouseButton.LEFT:
                yPicks = [self.dataUI.animationPicking.mousePositionInit[1],
                          self.dataUI.animationPicking.mousePosition[1]]
                xPicks = [self.dataUI.animationPicking.mousePositionInit[0],
                          self.dataUI.animationPicking.mousePosition[0]]
                if (self.dataUI.animationPicking.linePicking and
                        np.abs(yPicks[0]-yPicks[1]) > 1):
                    # The picking concerns multiple traces.
                    idMin, idMax = np.argmin(yPicks), np.argmax(yPicks)
                    minPick = np.max([0, np.ceil(yPicks[idMin])])
                    maxPick = np.min(
                        [np.floor(yPicks[idMax]), len(self.dataUI.sisData[0])-1])
                    for currPick in np.arange(minPick, maxPick+1, dtype=np.int16):
                        # Find x (time) for the given y (currPick).
                        m = (xPicks[idMax]-xPicks[idMin]) / \
                            (yPicks[idMax]-yPicks[idMin])
                        timePick = m*(currPick-yPicks[idMin]) + xPicks[idMin]
                        if timePick >= self.dataUI.beginTime[self.dataUI.sisFileId]:
                            self.ensure_manual_picking_mask()
                            self.dataUI.picking[self.dataUI.sisFileId,
                                                currPick] = timePick
                            # max(timePick*DEFAULT_ERROR, 0.000001) # Default error is 3%
                            self.dataUI.pickingError[self.dataUI.sisFileId,
                                                     currPick] = DEFAULT_ERROR
                            self.dataUI.manualPickingMask[
                                self.dataUI.sisFileId, currPick] = True
                elif not self.dataUI.animationPicking.linePicking:
                    self.select_picks_in_rectangle(xPicks, yPicks)
                self.dataUI.animationPicking.changedSelect = True
            # If right click and not dragging accross the pannel
            if event.button == MouseButton.RIGHT and ((time.time() - self.dataUI.animationPicking.timeOnClick) < self.dataUI.animationPicking.maxClickLength):
                if not (np.isnan(self.dataUI.picking[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect])):
                    if self.dataUI.animationPicking.mousePosition[0] > 0:
                        self.dataUI.pickingError[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect] = np.abs(
                            self.dataUI.animationPicking.mousePosition[0]-self.dataUI.picking[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect])/self.dataUI.picking[self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect]  # Default error is 3%
                        self.dataUI.animationPicking.changedSelect = True
        else:
            # Update the graph anyway for zoom in and pan updates
            self.dataUI.animationPicking.changedSelect = True
        # Picking is done for now, waiting for next click
        self.dataUI.animationPicking.notPicking = True
        self.dataUI.animationPicking.linePicking = False
        return 0

    def onKeyPress(self, event):
        """Handle Matplotlib key events reserved for future picker shortcuts."""
        print(event.key)
        # if event.key == "left" or event.key == "down":
        #     print('')

    def onPressModelling(self, event):
        """Select a movable hodograph control point near a mouse press."""
        if event.button == MouseButton.LEFT:
            sourceId = self.sourceSelector.currentIndex()
            receiversOrientation = self.receiversSelector.currentIndex()
            points = np.asarray(
                self.dataUI.modellingData.hodoPoints[sourceId][receiversOrientation])
            nbPoints = points.shape[0]
            mousePosition = np.repeat(np.asarray(
                [event.xdata, event.ydata]).reshape((1, 2)), nbPoints, axis=0)
            diff = np.abs(points-mousePosition)
            diffXBool = diff[:, 0] < self.dataUI.modellingAnimation.tol[0]
            diffYBool = diff[:, 1] < self.dataUI.modellingAnimation.tol[1]
            ptsId = np.logical_and(diffXBool, diffYBool)
            ptsId = np.where(ptsId)[0]
            if len(ptsId) > 1:
                ptsId = ptsId[-1]  # We take the last from both possibilities
            if len(ptsId) > 0:
                if ptsId > 0:
                    self.dataUI.modellingAnimation.currPointId = ptsId
                    self.dataUI.modellingAnimation.changingPts = True
                else:
                    self.dataUI.modellingAnimation.changingPts = False
        return

    def onReleaseModelling(self, event):
        """End the current hodograph control-point drag."""
        self.dataUI.modellingAnimation.changingPts = False

    def changeMouseModelling(self, event):  # For the modelling window
        """Move the selected hodograph point with the cursor during a drag."""
        if (event.inaxes is not None) and self.dataUI.modellingAnimation.changingPts:
            sourceId = self.sourceSelector.currentIndex()
            receiversOrientation = self.receiversSelector.currentIndex()
            ptsId = self.dataUI.modellingAnimation.currPointId
            self.dataUI.modellingData.hodoPoints[sourceId][receiversOrientation][ptsId, :] = [
                event.xdata, event.ydata]
        return 0

    def animationZoom(self, i):
        """Redraw the detail view around the cursor and current first-arrival pick."""
        axZoom = self.zoomGraph.axes
        # Get axis variables
        deltaT = float(
            self.dataUI.sisData[self.dataUI.sisFileId][0].stats.delta)
        nbPoints = self.dataUI.sisData[self.dataUI.sisFileId][0].stats.npts
        timeSEG2 = (self.dataUI.beginTime[self.dataUI.sisFileId] +
                    np.arange(nbPoints) * deltaT)
        # Change plot to go at the correct position:
        axZoom.clear()
        idx = np.greater_equal(timeSEG2, self.dataUI.animationPicking.mousePosition[0]-100*deltaT) & np.less_equal(
            timeSEG2, self.dataUI.animationPicking.mousePosition[0]+100*deltaT)
        timeZoom = timeSEG2[idx]
        axZoom.axhline(y=0, linewidth=0.5, color=[0.5, 0.5, 0.5])
        axZoom.axvline(x=0, linewidth=0.5, color=[0.5, 0.5, 0.5])
        vals = self.dataUI.sisData[self.dataUI.sisFileId][self.dataUI.animationPicking.currSelect].data[idx]
        trace_id = self.dataUI.animationPicking.currSelect
        isDefective = self.effective_defective_mask(
            self.dataUI.sisFileId)[trace_id]
        zeroOffsetArray = np.asarray(self.dataUI.zeroOffsetTraces)
        isZeroOffset = (zeroOffsetArray.ndim == 2 and
                        self.dataUI.sisFileId < zeroOffsetArray.shape[0] and
                        trace_id < zeroOffsetArray.shape[1] and
                        zeroOffsetArray[self.dataUI.sisFileId, trace_id])
        traceColor = '0.55' if isDefective else ('steelblue' if isZeroOffset else 'k')
        traceStyle = '--' if isDefective else (':' if isZeroOffset else '-')
        axZoom.plot(timeZoom, vals, color=traceColor, linestyle=traceStyle)
        if isDefective:
            reason = self.defective_trace_reason(
                self.dataUI.sisFileId, trace_id)
            label = 'DEFECTIVE GEOPHONE' + (f': {reason}' if reason else '')
            axZoom.text(0.02, 0.92, label, color='darkred',
                        transform=axZoom.transAxes, fontsize=9,
                        verticalalignment='top')
        elif isZeroOffset:
            axZoom.text(0.02, 0.92, 'ZERO-OFFSET TRACE (t0 reference)',
                        color='steelblue', transform=axZoom.transAxes,
                        fontsize=9, verticalalignment='top')
        axZoom.set_xlim(left=self.dataUI.animationPicking.mousePosition[0]-100 *
                        deltaT, right=self.dataUI.animationPicking.mousePosition[0]+100*deltaT)
        if len(vals) < 1:
            maxVal = 0.5
        else:
            maxVal = np.percentile(np.abs(vals), 75)
        axZoom.set_ylim(top=maxVal, bottom=-maxVal)  # autoscale(axis='y')
        # z = axZoom.get_ylim()
        axZoom.axvline(
            self.dataUI.animationPicking.mousePosition[0], color='r')
        currPicking = self.dataUI.picking[self.dataUI.sisFileId,
                                          self.dataUI.animationPicking.currSelect]
        if not (np.isnan(currPicking)):
            relativeError = self.dataUI.pickingError[
                self.dataUI.sisFileId, self.dataUI.animationPicking.currSelect]
            if not np.isfinite(relativeError):
                relativeError = DEFAULT_ERROR
            currError = relativeError * currPicking
            isSelected = (self.dataUI.animationPicking.currSelect in
                          self.dataUI.animationPicking.selectedPicks)
            pickColor = ('dodgerblue' if isSelected else
                         'darkorange' if relativeError > 0.05 else 'g')
            pickWidth = 2.5 if isSelected else 1.0
            axZoom.axvline(currPicking, color=pickColor, linewidth=pickWidth)
            axZoom.axvline(currPicking - currError, linestyle=':', color=pickColor)
            axZoom.axvline(currPicking + currError, linestyle=':', color=pickColor)
        axZoom.set_frame_on(False)
        axZoom.tick_params(
            axis='both',       # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            left=False,
            right=False,
            labelleft=False,
            labelbottom=False)  # labels along the bottom edge are off
        # draw_idle coalesces rapid mouse-motion redraw requests and keeps Qt
        # responsive while the user moves through the traces.
        self.zoomGraph.draw_idle()
        return 0

    def animationMain(self, i):
        """Redraw the visible portion of the active seismic gather when needed."""
        axMain = self.mainGraph.axes
        # Get axis variables
        deltaT = float(
            self.dataUI.sisData[self.dataUI.sisFileId][0].stats.delta)
        nbPoints = self.dataUI.sisData[self.dataUI.sisFileId][0].stats.npts
        timeSEG2 = (self.dataUI.beginTime[self.dataUI.sisFileId] +
                    np.arange(nbPoints) * deltaT)
        if not (self.dataUI.animationPicking.notPicking) or self.dataUI.animationPicking.changedSelect:
            # Change red graph + pick
            if not (self.dataUI.animationPicking.first):
                limitsY = axMain.get_ylim()
                limitsX = axMain.get_xlim()
            else:
                limitsX = self.mainGraph.homeXLimits

            # Only send samples inside the visible time window to Matplotlib.
            # This is substantially faster than plotting a complete long trace
            # and merely hiding most of it with set_xlim().
            visible = np.logical_and(
                timeSEG2 >= limitsX[0], timeSEG2 <= limitsX[1])
            if not np.any(visible):
                visible[np.argmin(np.abs(timeSEG2 - np.mean(limitsX)))] = True
            timeVisible = timeSEG2[visible]
            axMain.clear()
            i = 0
            for tr in self.dataUI.sisData[self.dataUI.sisFileId]:
                data = np.asarray(tr.data)[visible]
                amplitude = np.ptp(data)
                if not np.isfinite(amplitude) or amplitude == 0:
                    amplitude = 1.0
                data = (data/amplitude)*self.dataUI.visualAmplitudeScale+i
                isDefective = self.effective_defective_mask(
                    self.dataUI.sisFileId)[i]
                zeroOffsetArray = np.asarray(self.dataUI.zeroOffsetTraces)
                isZeroOffset = (zeroOffsetArray.ndim == 2 and
                                self.dataUI.sisFileId < zeroOffsetArray.shape[0] and
                                i < zeroOffsetArray.shape[1] and
                                zeroOffsetArray[self.dataUI.sisFileId, i])
                if isDefective:
                    axMain.plot(timeVisible, data, color='0.65', linestyle='--')
                elif isZeroOffset:
                    axMain.plot(timeVisible, data, color='steelblue', linestyle=':')
                elif i == self.dataUI.animationPicking.currSelect:
                    axMain.plot(timeVisible, data, color='r')
                else:
                    axMain.plot(timeVisible, data, color='k')
                axMain.axhline(y=i, linewidth=0.5, color=[0.5, 0.5, 0.5])
                currPicking = self.dataUI.picking[self.dataUI.sisFileId, i]
                if not (np.isnan(currPicking)):
                    relativeError = self.dataUI.pickingError[
                        self.dataUI.sisFileId, i]
                    if not np.isfinite(relativeError):
                        relativeError = DEFAULT_ERROR
                    currError = relativeError * currPicking
                    isSelected = i in self.dataUI.animationPicking.selectedPicks
                    pickColor = ('dodgerblue' if isSelected else
                                 'darkorange' if relativeError > 0.05 else 'g')
                    pickWidth = 2.5 if isSelected else 1.0
                    axMain.plot([currPicking, currPicking],
                                [i-0.5, i+0.5], color=pickColor,
                                linewidth=pickWidth)
                    axMain.plot(
                        [currPicking - currError, currPicking - currError],
                        [i-0.25, i+0.25], linestyle=':', color=pickColor)
                    axMain.plot(
                        [currPicking + currError, currPicking + currError],
                        [i-0.25, i+0.25], linestyle=':', color=pickColor)
                i += 1
            # axMain.axvline(x=0, linewidth=0.5, color=[0.5, 0.5, 0.5])
            axMain.xaxis.grid(True)
            axMain.set_xlim(left=limitsX[0], right=limitsX[1])
            if not (self.dataUI.animationPicking.first):
                axMain.set_ylim(bottom=limitsY[0], top=limitsY[1])
            else:  # First time plotting the graph, reset the home view
                self.mainGraph.homeYLimits = axMain.get_ylim()
            self.dataUI.animationPicking.first = False
            self.dataUI.animationPicking.changedSelect = False
            if not (self.dataUI.animationPicking.notPicking) and ((time.time() - self.dataUI.animationPicking.timeOnClick) > self.dataUI.animationPicking.maxClickLength):
                x = [self.dataUI.animationPicking.mousePositionInit[0],
                     self.dataUI.animationPicking.mousePosition[0]]
                y = [self.dataUI.animationPicking.mousePositionInit[1],
                     self.dataUI.animationPicking.mousePosition[1]]
                if self.dataUI.animationPicking.linePicking:
                    axMain.plot(x, y, color='g', linewidth=1.0)
                else:
                    selectionRectangle = Rectangle(
                        (min(x), min(y)), abs(x[1]-x[0]), abs(y[1]-y[0]),
                        edgecolor='dodgerblue', facecolor='dodgerblue',
                        alpha=0.18, linewidth=1.2)
                    axMain.add_patch(selectionRectangle)
            axMain.set_xlabel('Time [s]')
            axMain. set_ylabel('Trace #')
            # Do not redraw the full gather on every animation tick when
            # nothing changed; this was the main source of sluggishness.
            self.mainGraph.draw_idle()
        return 0

    # UI objects definition
    def _openGeometry(self):
        """Load a geometry file and the seismic files it references."""
        self.statusBar.showMessage('Openning Geometry file . . .')
        # Opening the geometry file:
        # The first argument returned is the filename and path
        fname, _ = QFileDialog.getOpenFileName(
            self, 'Open geometry file', filter='Geometry file (*.geometry)')
        if fname != "":
            headTail = os.path.split(fname)
            path = headTail[0]
            file = headTail[1]
            # nameSave = file[:-9] # remove the *.geometry
            # Retreive the datafiles names:
            SEG2Files = []
            SourcePosition = []
            ReceiversPosition = []
            sources = False
            receivers = False
            with open(fname) as f:
                Lines = f.read().splitlines()
            for line in Lines:
                if len(line.strip(' \t')) != 0:  # The line is stripped of spaces and tab
                    if line.startswith("SOURCES"):
                        sources = True
                        receivers = False
                    elif line.startswith("RECEIVERS"):
                        sources = False
                        receivers = True
                    else:
                        if sources:
                            # For robustness --> either tabs or spaces
                            tmp = re.split(r'  +|\t+', line.strip(' \t'))
                            name = tmp[0]
                            CurrSource = [float(i) for i in tmp[1:]]
                            SEG2Files.append(name)
                            SourcePosition.append(CurrSource)
                        elif receivers:
                            # For robustness --> either tabs or spaces
                            CurrReceiver = [float(i) for i in re.split(
                                r'  +|\t+', line.strip(' \t'))]
                            ReceiversPosition.append(CurrReceiver)
            # Check if Sources in List of Receivers --> Constitute Sensors array for output file:
            sensors = deepcopy(ReceiversPosition)
            sourcesId = np.zeros((len(SourcePosition),))
            i = 0
            for source in SourcePosition:
                # If the source is not in the receivers array
                if not (ReceiversPosition.count(source) == 1):
                    sensors.append(source)
                # We store the position of the Id position of the current source in the sensor array
                sourcesId[i] = sensors.index(source)
                i += 1

            if self.dataUI.dataLoaded == True:
                reply = QMessageBox.question(
                    self, 'Overwritting data . . .', 'Are you sure you want to overwrite the current dataset?', QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
                if reply == QMessageBox.Yes:
                    # Setting up the paths:
                    self.saveDataUI(path, file, SEG2Files,
                                    ReceiversPosition, sensors, sourcesId)

                    # Updating status bar
                    self.dataUI.dataLoaded = True
                    self.statusBar.showMessage(
                        f'{len(SEG2Files)} data files retreived from the geometry file with {len(sensors)} sensors.', 10000)
                else:
                    self.statusBar.showMessage(f'No data loaded', 2000)
            else:
                self.saveDataUI(path, file, SEG2Files,
                                ReceiversPosition, sensors, sourcesId)

                # Updating status bar
                self.dataUI.dataLoaded = True
                self.statusBar.showMessage(
                    f'{len(SEG2Files)} data files retreived from the geometry file with {len(sensors)} sensors.', 10000)

            # Return to the picking tab
            self.updateTab0()
            self.tabs.setCurrentIndex(0)
        else:
            self.statusBar.showMessage('No file loaded!', 2000)

    def saveDataUI(self, path, file, SEG2Files, ReceiversPosition, sensors, sourcesId):
        """Populate application state from parsed geometry and SEG-2/SEG-Y files."""
        # Setting up the paths:
        self.dataUI.paths.directory = path
        self.dataUI.paths.geometryFile = file
        self.dataUI.paths.nbFiles = len(sourcesId)
        self.dataUI.paths.seg2Files = SEG2Files

        # Setting up the data:
        self.dataUI.geometry.sensors = sensors
        self.dataUI.geometry.sourcesId = sourcesId
        self.dataUI.geometry.receivers = ReceiversPosition

        # Reading the datasets:
        for name in SEG2Files:
            _, ext = os.path.splitext(name)
            if ext == '.segy' or ext == '.sgy':
                st = read(os.path.join(path, name), 'SEGY')
                self.dataUI.beginTime.append(0)
            elif ext == '.seg2' or ext == '.sg2':
                st = read(os.path.join(path, name), 'SEG2')
                self.dataUI.beginTime.append(
                    float(st[0].stats.seg2["DELAY"].replace(',', '.')))
            else:
                st = read(os.path.join(path, name), 'SEG2')
                self.dataUI.beginTime.append(
                    float(st[0].stats.seg2["DELAY"].replace(',', '.')))
            # If the number of geophones does not match between the loaded array and the gemoetry
            if len(st) != len(ReceiversPosition):
                raise Exception(
                    'The file referenced in the geometry file does not match the geometry of the array!')
            self.dataUI.sisData.append(st)

        self.dataUI.sisDataOriginal = deepcopy(self.dataUI.sisData)
        # self.filterSignal()
        self.dataUI.picking = np.empty(
            (len(self.dataUI.paths.seg2Files), len(self.dataUI.sisData[0])))
        self.dataUI.picking[:] = np.nan
        self.dataUI.pickingError = np.empty(
            (len(self.dataUI.paths.seg2Files), len(self.dataUI.sisData[0])))
        self.dataUI.pickingError[:] = np.nan
        self.dataUI.manualPickingMask = np.zeros_like(
            self.dataUI.picking, dtype=bool)
        self.dataUI.defectiveGeophones = np.zeros(
            len(self.dataUI.sisData[0]), dtype=bool)
        self.dataUI.defectiveGeophoneReasons = [
            '' for _ in range(len(self.dataUI.sisData[0]))]
        self.dataUI.manualDefectiveOverrides = np.zeros_like(
            self.dataUI.picking, dtype=np.int8)
        self.dataUI.zeroOffsetTraces = np.zeros_like(
            self.dataUI.picking, dtype=bool)
        self.spinBoxCurrSelect.setMinimum(0)
        self.spinBoxCurrSelect.setMaximum(len(self.dataUI.sisData[0])-1)
        self.spinBoxCurrSelect.setPrefix('Trace number ')
        self.spinBoxCurrSelect.valueChanged.connect(self.traceNumberChanged)
        self.spinBoxCurrSelect.setEnabled(True)
        self.buttonManageDefective.setEnabled(True)

    def fftGraph(self):
        """Toggle the main gather between time-domain traces and FFT spectra."""
        if self.dataUI.dataLoaded:
            if not (self.dataUI.animationPicking.fftShowed):
                axMain = self.mainGraph.axes
                # We stop the animation:
                self.aniZoom.event_source.stop()
                self.aniMain.event_source.stop()
                self.mainGraph.mpl_disconnect(self.connectMouse)
                self.mainGraph.mpl_disconnect(self.connectPress)
                self.mainGraph.mpl_disconnect(self.connectRelease)
                self.mainGraph.mpl_disconnect(self.connectKeyPress)
                self.spinBoxCurrSelect.setEnabled(False)
                self.buttonManageDefective.setEnabled(False)
                # Changing the mainGraph to show fft instead of traces.
                axMain.clear()
                deltaT = float(
                    self.dataUI.sisData[self.dataUI.sisFileId][0].stats.delta)
                nbPoints = self.dataUI.sisData[self.dataUI.sisFileId][0].stats.npts
                freq = np.fft.fftfreq(nbPoints, deltaT)
                for i, tr in enumerate(self.dataUI.sisData[self.dataUI.sisFileId]):
                    data = tr.data
                    fftData = np.abs(np.fft.fft(data))
                    fftData = fftData/max(fftData) + i
                    axMain.axhline(y=i, linewidth=0.5, color=[0.5, 0.5, 0.5])
                    if i == self.dataUI.animationPicking.currSelect:
                        axMain.plot(freq[:nbPoints//4],
                                    fftData[:nbPoints//4], color='r')
                        axMain.fill_between(
                            freq[:nbPoints//4], np.ones_like(freq[:nbPoints//4])*i, fftData[:nbPoints//4], color='r')
                    else:
                        axMain.plot(freq[:nbPoints//4],
                                    fftData[:nbPoints//4], color='k')
                        axMain.fill_between(
                            freq[:nbPoints//4], np.ones_like(freq[:nbPoints//4])*i, fftData[:nbPoints//4], color='k')
                axMain.xaxis.grid(True)
                axMain.set_xlabel('Frequency [Hz]')
                axMain. set_ylabel('Trace #')
                self.mainGraph.draw()
                self.dataUI.animationPicking.fftShowed = True
                # Change the menu item to get the original graph back
                fftAction = self.dataUI.animationPicking.fftAnim
                fftAction.setStatusTip(
                    'Return to the original time-series graphs')
                fftAction.setText('Show the time-series dataset')
                self.mainGraph.homeXLimits = axMain.get_xlim()
                self.mainGraph.homeYLimits = axMain.get_ylim()
            else:
                # We begin back the animation:
                self.aniZoom.event_source.start()
                self.aniMain.event_source.start()
                self.connectMouse = self.mainGraph.mpl_connect(
                    'motion_notify_event', self.changeMouse)
                self.connectPress = self.mainGraph.mpl_connect(
                    'button_press_event', self.onPress)
                self.connectRelease = self.mainGraph.mpl_connect(
                    'button_release_event', self.onRelease)
                self.connectKeyPress = self.mainGraph.mpl_connect(
                    'key_press_event', self.onKeyPress)
                self.spinBoxCurrSelect.setEnabled(True)
                self.buttonManageDefective.setEnabled(True)
                self.dataUI.animationPicking.first = True
                self.dataUI.animationPicking.changedSelect = True
                self.dataUI.animationPicking.fftShowed = False
                # Change the menu item to get the original graph back
                fftAction = self.dataUI.animationPicking.fftAnim
                fftAction.setStatusTip(
                    'Show the FFT transform of the different traces')
                fftAction.setText('Show FFT of dataset')

    def dcFilter(self):
        """Remove the constant (DC) component from every loaded trace."""
        for st in self.dataUI.sisData:
            for tr in st:
                tr.detrend('constant')
        self.dataUI.animationPicking.changedSelect = True

    def highpassFilter(self):
        """Prompt for and apply a high-pass filter to every loaded trace."""
        band, ok = QInputDialog.getDouble(
            self, "High-pass filter", "Frequency [Hz]", 5.0, 0.0, 10000.0)
        if ok:
            for st in self.dataUI.sisData:
                for tr in st:
                    tr.filter('highpass', freq=band)
        self.dataUI.animationPicking.changedSelect = True

    def lowpassFilter(self):
        """Prompt for and apply a low-pass filter to every loaded trace."""
        band, ok = QInputDialog.getDouble(
            self, "Low-pass filter", "Frequency [Hz]", 100.0, 0.0, 10000.0)
        if ok:
            for st in self.dataUI.sisData:
                for tr in st:
                    tr.filter('lowpass', freq=band)
        self.dataUI.animationPicking.changedSelect = True

    def trimFilter(self):
        """Prompt for a time interval and trim all loaded streams to it."""
        deltaT = float(
            self.dataUI.sisData[self.dataUI.sisFileId][0].stats.delta)
        nbPoints = self.dataUI.sisData[self.dataUI.sisFileId][0].stats.npts
        tEnd = self.dataUI.beginTime[self.dataUI.sisFileId]+nbPoints*deltaT
        trim, ok = QInputDialog.getDouble(
            self, "Trimming", "Time: ", tEnd/2, 0.0, tEnd, 2)
        if ok:
            for st in self.dataUI.sisData:
                st.trim(starttime=st[0].stats.starttime,
                        endtime=st[0].stats.starttime + trim)
        self.dataUI.animationPicking.changedSelect = True

    def resetFilters(self):
        """Restore the unfiltered seismic streams saved at load time."""
        self.dataUI.sisData = deepcopy(self.dataUI.sisDataOriginal)
        self.dataUI.animationPicking.changedSelect = True

    def filterSignal(self, dcFilter: bool = True, highpass: float = 5.0, lowpass: float = 100.0):
        """Apply the standard detrend, high-pass, and low-pass preprocessing."""
        self.dataUI.sisData = deepcopy(self.dataUI.sisDataOriginal)
        if not (dcFilter) and (highpass > 0.0 or lowpass < np.Inf):
            QMessageBox.warning(
                self, 'Warning !', 'Using the lowpass and/or highpass filter\n without the DC filter is not advised!')
        for st in self.dataUI.sisData:
            for tr in st:
                if dcFilter:
                    tr.detrend('constant')
                if highpass > 0.0:
                    tr.filter('highpass', freq=highpass)
                if lowpass < np.Inf:
                    tr.filter('lowpass', freq=lowpass)
        self.dataUI.animationPicking.changedSelect = True

    def detect_defective_geophones(self):
        """Detect persistently dead, clipped, invalid, or noisy channels.

        Returns
        -------
        defective : numpy.ndarray of bool
            Per-trace mask. A trace is marked when anomalous behaviour appears
            in at least half of the loaded shots.
        reasons : list of str
            Human-readable diagnostic reasons for every trace.

        Notes
        -----
        This method reads all streams in ``self.dataUI.sisData`` but does not
        modify them. It combines amplitude, flatness, clipping, impulsiveness,
        and high-frequency roughness tests across shots.
        """
        if len(self.dataUI.sisData) == 0:
            return np.array([], dtype=bool), []
        nb_files = len(self.dataUI.sisData)
        nb_traces = len(self.dataUI.sisData[0])
        votes = np.zeros(nb_traces, dtype=int)
        reason_votes = [dict() for _ in range(nb_traces)]

        for stream in self.dataUI.sisData:
            scales = np.full(nb_traces, np.nan)
            roughnesses = np.full(nb_traces, np.nan)
            preliminary = [[] for _ in range(nb_traces)]
            for trace_id, trace in enumerate(stream):
                values = np.asarray(trace.data, dtype=float)
                finite = np.isfinite(values)
                if values.size < 24 or np.mean(finite) < 0.99:
                    preliminary[trace_id].append('invalid samples')
                    continue
                values = values[finite]
                centred = values-np.median(values)
                scale = np.median(np.abs(centred))/0.67448975
                if not np.isfinite(scale) or scale <= np.finfo(float).eps:
                    preliminary[trace_id].append('dead/constant')
                    continue
                scales[trace_id] = scale

                differences = np.diff(values)
                roughnesses[trace_id] = (
                    np.median(np.abs(differences)) /
                    max(np.median(np.abs(centred)), np.finfo(float).eps))
                flat_fraction = (np.mean(differences == 0.0)
                                 if differences.size else 1.0)
                minimum = np.min(values)
                maximum = np.max(values)
                clipped_fraction = max(
                    np.mean(values == minimum), np.mean(values == maximum))
                percentile_99 = np.percentile(np.abs(centred), 99)
                crest = np.max(np.abs(centred))/max(
                    percentile_99, np.finfo(float).eps)
                if flat_fraction > 0.995:
                    preliminary[trace_id].append('mostly constant')
                if clipped_fraction > 0.03:
                    preliminary[trace_id].append('clipped/saturated')
                if crest > 30.0:
                    preliminary[trace_id].append('impulsive spikes')

            valid_scales = scales[np.isfinite(scales) & (scales > 0)]
            reference_scale = (float(np.median(valid_scales))
                               if valid_scales.size else np.nan)
            valid_roughness = roughnesses[np.isfinite(roughnesses)]
            reference_roughness = (float(np.median(valid_roughness))
                                   if valid_roughness.size else np.nan)
            for trace_id in range(nb_traces):
                reasons = preliminary[trace_id]
                if np.isfinite(reference_scale) and np.isfinite(scales[trace_id]):
                    ratio = scales[trace_id]/reference_scale
                    if ratio < 0.002:
                        reasons.append('abnormally weak')
                    elif ratio > 50.0:
                        reasons.append('abnormally noisy')
                if (np.isfinite(reference_roughness) and
                        np.isfinite(roughnesses[trace_id]) and
                        roughnesses[trace_id] > max(
                            0.20, 3.0*reference_roughness)):
                    reasons.append('incoherent high-frequency noise')
                if reasons:
                    votes[trace_id] += 1
                    for reason in set(reasons):
                        reason_votes[trace_id][reason] = (
                            reason_votes[trace_id].get(reason, 0)+1)

        # A hardware problem should normally persist over several shots. For a
        # single-file dataset the conservative per-gather tests are used.
        minimum_votes = max(1, int(np.ceil(0.5*nb_files)))
        defective = votes >= minimum_votes
        reasons = []
        for trace_id in range(nb_traces):
            ordered = sorted(
                reason_votes[trace_id].items(), key=lambda item: item[1], reverse=True)
            reasons.append(', '.join(reason for reason, _ in ordered))
        return defective, reasons

    def remove_spatial_outliers(self, selected, candidates_by_trace,
                                branches, anchored, excluded,
                                receiver_distances, dt):
        """Reject or replace picks outside the robust segmented trajectory.

        Parameters
        ----------
        selected : dict
            Current trace-to-candidate selection; it is updated in place.
        candidates_by_trace : sequence or mapping
            All waveform candidates indexed by trace identifier.
        branches : sequence of array-like
            Trace identifiers grouped and ordered by profile side.
        anchored, excluded : array-like of bool
            Masks for protected manual anchors and unusable traces.
        receiver_distances : array-like
            Source-receiver distances in metres.
        dt : float
            Sampling interval in seconds.

        Returns
        -------
        dict
            Updated selection. Large-residual picks are replaced by a credible
            nearby alternative or removed when no such alternative exists.
        """
        min_tolerance = max(12.0, 0.002/max(dt, np.finfo(float).eps))
        for branch in branches:
            usable = [int(i) for i in branch
                      if not excluded[int(i)] and int(i) in selected]
            for trace_id in usable:
                if anchored[trace_id]:
                    continue
                candidate = selected[trace_id]
                prediction = candidate.get('spatial_prediction')
                if prediction is None:
                    continue
                residual = abs(candidate['index']-prediction)
                if residual <= 2.5*min_tolerance:
                    continue
                alternatives = [
                    item for item in candidates_by_trace[trace_id]
                    if abs(item['index']-prediction) <= 2.0*min_tolerance and
                    item['confidence'] >= 0.30]
                if alternatives:
                    replacement = min(
                        alternatives,
                        key=lambda item: item['cost'] +
                        0.45*abs(item['index']-prediction)/min_tolerance)
                    replacement = dict(replacement)
                    replacement['spatial_prediction'] = float(prediction)
                    replacement['spatial_residual'] = float(
                        abs(replacement['index']-prediction))
                    replacement['spatial_agreement'] = float(np.exp(
                        -replacement['spatial_residual']/min_tolerance))
                    replacement['segment_id'] = candidate.get('segment_id', 0)
                    selected[trace_id] = replacement
                else:
                    selected.pop(trace_id, None)
        return selected

    def validate_branch_anchor(self, candidates_by_trace, branch,
                               anchored, excluded, receiver_distances, dt):
        """Validate and lock the closest usable receiver before propagation.

        Parameters
        ----------
        candidates_by_trace : sequence or mapping
            Candidate dictionaries indexed by trace identifier. A validated
            candidate is locked by assigning it a negative cost.
        branch : array-like of int
            Trace identifiers ordered away from the source.
        anchored, excluded : array-like of bool
            Masks for existing trusted picks and unusable traces.
        receiver_distances : array-like
            Source-receiver distances in metres.
        dt : float
            Sampling interval in seconds.

        Returns
        -------
        valid : bool
            Whether a sufficiently strong and coherent anchor was found.
        trace_id : int or None
            Identifier of the validated anchor, or the trace requiring manual
            review when validation fails.
        """
        usable = [int(i) for i in branch if not excluded[int(i)]]
        if not usable:
            return False, None
        nearest = usable[0]
        if anchored[nearest]:
            return True, nearest
        if not candidates_by_trace[nearest]:
            return False, nearest

        # Determine the origin from a short coherent sequence rather than from
        # the strongest isolated candidate on the nearest trace.
        provisional = select_consistent_candidates(
            candidates_by_trace, usable[:min(6, len(usable))],
            positions=receiver_distances,
            base_tolerance=max(8.0, 0.0015/max(dt, np.finfo(float).eps)),
            regularize_piecewise=False)
        if nearest not in provisional:
            return False, nearest
        anchor = provisional[nearest]
        strong_enough = (anchor['confidence'] >= 0.50 and
                         anchor.get('snr', 0.0) >= 3.0)

        positive_distances = receiver_distances[receiver_distances > 0]
        zero_tolerance = (max(1e-9, 1e-6*np.min(positive_distances))
                          if positive_distances.size else 1e-9)
        nearest_is_zero_offset = (
            receiver_distances[nearest] <= zero_tolerance)
        if nearest_is_zero_offset:
            # Source coupling and near-field energy make the zero-offset onset
            # fundamentally different from the first sloping travel-time
            # branch. Validate it from its own waveform, and validate the
            # following receivers separately instead of extrapolating their
            # line back through the source trace.
            anchor = max(
                candidates_by_trace[nearest],
                key=lambda item: item.get('confidence', 0.0) +
                0.05*np.clip(np.log(max(item.get('snr', 1.0), 1.0)) /
                             np.log(50.0), 0.0, 1.0))
            follower_ids = usable[1:min(7, len(usable))]
            follower_selection = select_consistent_candidates(
                candidates_by_trace, follower_ids,
                positions=receiver_distances,
                base_tolerance=max(
                    8.0, 0.0015/max(dt, np.finfo(float).eps)),
                regularize_piecewise=False)
            follower_ids = [i for i in follower_ids
                            if i in follower_selection]
            strong_follower_count = sum(
                follower_selection[i].get('confidence', 0.0) >= 0.60 and
                follower_selection[i].get('snr', 0.0) >= 3.0
                for i in follower_ids)
            followers_coherent = len(follower_ids) >= 3
            if followers_coherent:
                follower_x = np.asarray(
                    [receiver_distances[i] for i in follower_ids], dtype=float)
                follower_y = np.asarray(
                    [follower_selection[i]['index'] for i in follower_ids],
                    dtype=float)
                follower_prediction = _robust_line_prediction(
                    follower_x, follower_y)
                follower_residual = np.abs(
                    follower_y-follower_prediction)
                follower_tolerance = max(
                    24.0, 0.003/max(dt, np.finfo(float).eps))
                followers_coherent = (
                    np.median(follower_residual) <= follower_tolerance and
                    np.max(follower_residual) <= 2.5*follower_tolerance)

            zero_signal_strong = (
                anchor['confidence'] >= 0.55 and
                anchor.get('snr', 0.0) >= 5.0)
            zero_signal_very_strong = (
                anchor['confidence'] >= 0.75 and
                anchor.get('snr', 0.0) >= 10.0)
            neighbourhood_valid = (
                followers_coherent or strong_follower_count >= 3)
            if (zero_signal_strong and
                    (neighbourhood_valid or zero_signal_very_strong)):
                locked_anchor = dict(anchor)
                locked_anchor['cost'] = -25.0
                # Keep the t0 observation, but do not make the near-field
                # source pulse define the slope followed by other receivers.
                locked_anchor['exclude_from_spatial_trend'] = True
                candidates_by_trace[nearest] = [locked_anchor]
                return True, nearest

            # A doubtful source trace must not discard an otherwise clean
            # branch. Start at the first non-zero receiver and leave the
            # zero-offset trace unpicked for manual review.
            if followers_coherent:
                fallback_id = follower_ids[0]
                fallback = follower_selection[fallback_id]
                fallback_strong = (
                    fallback['confidence'] >= 0.50 and
                    fallback.get('snr', 0.0) >= 3.0)
                if fallback_strong:
                    locked_fallback = dict(fallback)
                    locked_fallback['cost'] = -25.0
                    candidates_by_trace[fallback_id] = [locked_fallback]
                    candidates_by_trace[nearest] = []
                    return True, fallback_id
            return False, nearest

        # The first receiver away from the source can legitimately belong to a
        # very steep near-field segment.  In that case it must not be rejected
        # merely because the next receivers already follow another slope.  The
        # dedicated near-source detector is deliberately allowed to override
        # this extrapolation only for an exceptionally clear, persistent onset;
        # this remains much stricter than accepting the best generic candidate.
        near_source_anchor = (
            anchor.get('near_source_onset', False) and
            anchor['confidence'] >= 0.78 and
            anchor.get('snr', 0.0) >= 12.0)
        if near_source_anchor:
            locked_anchor = dict(anchor)
            locked_anchor['cost'] = -25.0
            locked_anchor['near_source_anchor'] = True
            candidates_by_trace[nearest] = [locked_anchor]
            return True, nearest

        provisional_indices = [
            provisional[i]['index'] for i in usable[:min(6, len(usable))]
            if i in provisional]
        coherent = True
        if len(provisional_indices) >= 3:
            predicted_anchor = (2.0*provisional_indices[1] -
                                provisional_indices[2])
            local_step = abs(provisional_indices[2]-provisional_indices[1])
            coherent = abs(anchor['index']-predicted_anchor) <= max(
                16.0, 2.5*max(local_step, 1.0))
        if not (strong_enough and coherent):
            return False, nearest

        # Lock the validated onset so lateral regularisation cannot replace it
        # with a later, smoother phase.
        locked_anchor = dict(anchor)
        locked_anchor['cost'] = -25.0
        candidates_by_trace[nearest] = [locked_anchor]
        return True, nearest

    def refine_manual_candidate(self, trace, pick_idx, signal_start,
                                signal_stop, radius=100):
        """Refine a manual pick locally while retaining trusted-anchor status.

        Parameters
        ----------
        trace : obspy.Trace
            Seismic trace containing the manually picked arrival.
        pick_idx : int
            Current pick position as a sample index.
        signal_start, signal_stop : int
            Half-open valid signal interval in sample indices.
        radius : int, default=100
            Maximum search distance on either side of the manual pick.

        Returns
        -------
        dict
            Refined candidate metadata. Confidence is forced to one and cost
            to a negative value so later spatial regularization preserves it.
        """
        local_start = max(signal_start, int(pick_idx)-radius)
        local_stop = min(signal_stop, int(pick_idx)+radius+1)
        candidates = hybrid_pick_candidates(
            trace.data, local_start, local_stop, max_candidates=6)
        if not candidates:
            return {
                'index': int(pick_idx), 'confidence': 1.0,
                'cost': -50.0, 'snr': np.inf}
        refined = min(
            candidates,
            key=lambda item: item['cost'] +
            0.35*abs(item['index']-pick_idx)/max(radius, 1))
        refined = dict(refined)
        refined['confidence'] = 1.0
        refined['cost'] = -50.0
        return refined

    def improve_zero_offset_candidates(self, candidates_by_trace, stream,
                                       zero_offset, defective, anchored, distances,
                                       signal_start, signal_stop):
        """Guide zero-offset picks from nearby arrivals and refine them locally.

        Parameters
        ----------
        candidates_by_trace : sequence or mapping
            Mutable candidate lists indexed by trace identifier.
        stream : obspy.Stream
            Seismic traces for the current shot.
        zero_offset, defective, anchored : array-like of bool
            Masks identifying colocated traces, unusable channels, and trusted
            manual picks.
        distances : array-like
            Source-receiver distances in metres.
        signal_start, signal_stop : int
            Half-open valid signal interval in sample indices.

        Returns
        -------
        None

        Notes
        -----
        The function extrapolates nearby candidate indices to zero distance and
        promotes the best local waveform candidate in each zero-offset trace.
        ``candidates_by_trace`` is modified in place.
        """
        zero_ids = np.flatnonzero(zero_offset & ~defective)
        neighbour_ids = np.flatnonzero(~zero_offset & ~defective)
        neighbour_ids = neighbour_ids[np.argsort(distances[neighbour_ids])]
        for zero_id in zero_ids:
            if anchored[int(zero_id)]:
                continue
            used_distances = []
            used_indices = []
            for trace_id in neighbour_ids:
                candidates = candidates_by_trace[int(trace_id)]
                if not candidates:
                    continue
                candidate = max(
                    candidates,
                    key=lambda item: item['confidence'] +
                    0.10*np.clip(np.log(max(item.get('snr', 1.0), 1.0)) /
                                 np.log(50.0), 0.0, 1.0))
                if candidate['confidence'] >= 0.45:
                    used_distances.append(distances[trace_id])
                    used_indices.append(candidate['index'])
                if len(used_indices) >= 5:
                    break
            if len(used_indices) < 2 or len(np.unique(used_distances)) < 2:
                continue
            if len(used_indices) >= 3:
                theil_result = stats.theilslopes(
                    used_indices, used_distances)
                # Older SciPy releases return a plain tuple; newer releases
                # expose the same value through the ``intercept`` attribute.
                estimate = (theil_result.intercept
                            if hasattr(theil_result, 'intercept')
                            else theil_result[1])
            else:
                _, estimate = np.polyfit(used_distances, used_indices, 1)
            estimate = int(round(estimate))
            radius = 120
            local_candidates = hybrid_pick_candidates(
                stream[int(zero_id)].data,
                max(signal_start, estimate-radius),
                min(signal_stop, estimate+radius+1),
                max_candidates=6)
            if not local_candidates:
                continue
            refined = min(
                local_candidates,
                key=lambda item: item['cost'] +
                0.40*abs(item['index']-estimate)/radius)
            candidates_by_trace[int(zero_id)] = [refined] + [
                item for item in candidates_by_trace[int(zero_id)]
                if item['index'] != refined['index']]

    def autoPicking(self):
        """Automatically pick first arrivals with a hybrid multi-trace method.

        Returns
        -------
        None

        Notes
        -----
        The method operates on every loaded shot. It combines STA/LTA and AIC
        waveform candidates with receiver-to-receiver spatial consistency,
        excludes defective channels, treats zero-offset traces separately, and
        optionally preserves or refines existing manual picks. Accepted picks,
        relative errors, masks, and display state are written to ``self.dataUI``.
        A dialog requests a manual anchor when the first reliable trace of a
        profile branch cannot be validated automatically.
        """
        if not AIC_SIMPLE_IMPORT:
            QMessageBox.warning(
                self, 'Warning !', 'Auto-picking not implemented in current version of obspy.')
            return
        if not self.dataUI.dataLoaded or len(self.dataUI.sisData) == 0:
            QMessageBox.warning(self, 'Warning!', 'No seismic data are loaded.')
            return
        self.dataUI.animationPicking.selectedPicks.clear()

        self.ensure_manual_picking_mask()
        manual_existing = np.logical_and(
            np.isfinite(self.dataUI.picking), self.dataUI.manualPickingMask)
        manual_mode = 'ignore'
        if manual_existing.any():
            choices = [
                'Keep manual picks unchanged and use them as anchors',
                'Refine manual picks locally and use them as anchors',
                'Ignore manual picks and recompute everything']
            choice, ok = QInputDialog.getItem(
                self, 'Manual picks detected',
                'How should existing manual picks be handled?',
                choices, 0, False)
            if not ok:
                return
            manual_mode = ('keep' if choice == choices[0]
                           else 'refine' if choice == choices[1]
                           else 'ignore')

        requested_start = self.timeWindowStart.value()
        requested_end = self.timeWindowEnd.value()
        defective, defect_reasons = self.detect_defective_geophones()
        self.dataUI.defectiveGeophones = defective
        self.dataUI.defectiveGeophoneReasons = defect_reasons
        self.dataUI.zeroOffsetTraces = np.zeros_like(
            self.dataUI.picking, dtype=bool)
        # Count the effective mask below: it includes shot-specific user
        # decisions in addition to the automatic global diagnosis.
        defective_count = 0
        nb_traces = sum(len(stream) for stream in self.dataUI.sisData)
        progress = QProgressBar(self)
        progress.setRange(0, nb_traces)
        progress.setValue(0)
        self.statusBar.addWidget(progress, 1)
        self.statusBar.showMessage('Hybrid automatic picking in progress...')

        accepted = 0
        rejected = 0
        processed = 0
        zero_offset_count = 0
        manual_anchor_requests = []
        try:
            for file_id, stream in enumerate(self.dataUI.sisData):
                dt = float(stream[0].stats.delta)
                begin = float(self.dataUI.beginTime[file_id])
                npts = int(stream[0].stats.npts)
                signal_end = begin + (npts-1)*dt
                window_start = max(begin, requested_start)
                window_end = min(signal_end, requested_end)
                if window_end-window_start < 24*dt:
                    window_start, window_end = begin, signal_end
                start_idx = max(0, int(np.ceil((window_start-begin)/dt)))
                stop_idx = min(npts, int(np.floor((window_end-begin)/dt))+1)

                receivers = np.asarray(self.dataUI.geometry.receivers, dtype=float)
                source_id = int(self.dataUI.geometry.sourcesId[file_id])
                source = np.asarray(
                    self.dataUI.geometry.sensors[source_id], dtype=float)
                distances = calculateDistance(receivers, source)
                receiver_steps = np.linalg.norm(
                    np.diff(receivers, axis=0), axis=1)
                receiver_steps = receiver_steps[
                    np.isfinite(receiver_steps) & (receiver_steps > 0)]
                near_source_radius = (2.1*float(np.median(receiver_steps))
                                      if receiver_steps.size else 0.0)
                positive_distances = distances[distances > 0]
                zero_tolerance = (max(1e-9, 1e-6*np.min(positive_distances))
                                  if positive_distances.size else 1e-9)
                zero_offset = distances <= zero_tolerance
                self.dataUI.zeroOffsetTraces[file_id, :] = zero_offset
                # A zero-offset trace is a useful t0 reference. It may contain
                # a negative arrival before correction and must still be picked.
                excluded = self.effective_defective_mask(file_id)
                defective_count += int(np.sum(excluded))
                zero_offset_count += int(np.sum(zero_offset))

                candidates_by_trace = []
                anchored = np.zeros(len(stream), dtype=bool)
                for trace_id, trace in enumerate(stream):
                    old_pick = self.dataUI.picking[file_id, trace_id]
                    is_manual = bool(manual_existing[file_id, trace_id])
                    if excluded[trace_id]:
                        candidates = []
                    elif manual_mode == 'keep' and is_manual:
                        pick_idx = int(round((old_pick-begin)/dt))
                        candidates = [{
                            'index': int(np.clip(pick_idx, 0, npts-1)),
                            'confidence': 1.0,
                            'cost': -50.0,
                            'snr': np.inf}]
                        anchored[trace_id] = True
                    elif manual_mode == 'refine' and is_manual:
                        pick_idx = int(round((old_pick-begin)/dt))
                        candidates = [self.refine_manual_candidate(
                            trace, pick_idx, start_idx, stop_idx)]
                        anchored[trace_id] = True
                    else:
                        candidates = hybrid_pick_candidates(
                            trace.data, start_idx, stop_idx)
                        if distances[trace_id] <= near_source_radius:
                            onset = near_source_energy_candidate(
                                trace.data, begin, dt, start_idx, stop_idx)
                            if onset is not None:
                                by_index = {
                                    candidate['index']: candidate
                                    for candidate in candidates}
                                previous = by_index.get(onset['index'])
                                if (previous is None or
                                        onset['confidence'] >
                                        previous['confidence']):
                                    by_index[onset['index']] = onset
                                candidates = sorted(
                                    by_index.values(),
                                    key=lambda item: item['confidence'],
                                    reverse=True)[:7]
                    candidates_by_trace.append(candidates)
                    processed += 1
                    progress.setValue(processed)
                    if processed % 8 == 0:
                        QApplication.processEvents()

                self.improve_zero_offset_candidates(
                    candidates_by_trace, stream, zero_offset, excluded,
                    anchored,
                    distances, start_idx, stop_idx)
                branches = self.trace_branches(file_id)
                valid_branches = []
                for branch in branches:
                    anchor_valid, nearest_trace = self.validate_branch_anchor(
                        candidates_by_trace, branch, anchored, excluded,
                        distances, dt)
                    if anchor_valid:
                        valid_branches.append(branch)
                    elif nearest_trace is not None:
                        manual_anchor_requests.append(
                            (file_id, nearest_trace, window_start, window_end))
                selected = {}
                spatial_tolerance = max(
                    8.0, 0.0015/max(dt, np.finfo(float).eps))
                for branch in valid_branches:
                    selected.update(select_consistent_candidates(
                        candidates_by_trace, branch,
                        positions=distances,
                        base_tolerance=spatial_tolerance))
                selected = self.remove_spatial_outliers(
                    selected, candidates_by_trace, valid_branches,
                    anchored, excluded, distances, dt)

                # Combine the waveform score with agreement to the confirmed
                # piecewise-linear trajectory (including at slope breaks).
                for branch in valid_branches:
                    branch_selected = [i for i in branch if i in selected]
                    for trace_id in branch_selected:
                        candidate = selected[trace_id]
                        if anchored[trace_id]:
                            continue
                        agreement = candidate.get('spatial_agreement', 1.0)
                        candidate['confidence'] = float(
                            0.45*candidate['confidence'] + 0.55*agreement)

                for trace_id in range(len(stream)):
                    if excluded[trace_id]:
                        if not (manual_existing[file_id, trace_id] and
                                manual_mode in ('keep', 'refine')):
                            self.dataUI.picking[file_id, trace_id] = np.nan
                            self.dataUI.pickingError[file_id, trace_id] = np.nan
                            self.dataUI.manualPickingMask[file_id, trace_id] = False
                        continue
                    if (trace_id not in selected and
                            manual_existing[file_id, trace_id] and
                            manual_mode in ('keep', 'refine')):
                        # Never erase a user pick merely because automatic
                        # validation rejected the surrounding branch.
                        continue
                    if anchored[trace_id]:
                        if manual_mode == 'refine' and trace_id in selected:
                            refined = selected[trace_id]
                            self.dataUI.picking[file_id, trace_id] = (
                                begin+refined['index']*dt)
                            self.dataUI.pickingError[file_id, trace_id] = DEFAULT_ERROR
                            self.dataUI.manualPickingMask[file_id, trace_id] = True
                        accepted += 1
                        continue
                    candidate = selected.get(trace_id)
                    if candidate is None or candidate['confidence'] < 0.35:
                        self.dataUI.picking[file_id, trace_id] = np.nan
                        self.dataUI.pickingError[file_id, trace_id] = np.nan
                        self.dataUI.manualPickingMask[file_id, trace_id] = False
                        rejected += 1
                        continue
                    picked_time = begin + candidate['index']*dt
                    confidence = candidate['confidence']
                    relative_error = DEFAULT_ERROR*(1.0 + 2.0*(1.0-confidence))
                    self.dataUI.picking[file_id, trace_id] = picked_time
                    self.dataUI.pickingError[file_id, trace_id] = relative_error
                    self.dataUI.manualPickingMask[file_id, trace_id] = bool(
                        anchored[trace_id])
                    accepted += 1
        finally:
            self.statusBar.removeWidget(progress)
            progress.deleteLater()

        self.dataUI.pickingDone = accepted > 0
        self.dataUI.animationPicking.changedSelect = True
        if manual_anchor_requests:
            file_id, trace_id, window_start, window_end = manual_anchor_requests[0]
            file_name = os.path.basename(self.dataUI.paths.seg2Files[file_id])
            self.comboBoxFilesPicking.setCurrentIndex(file_id)
            self.spinBoxCurrSelect.setValue(trace_id)
            self.dataUI.animationPicking.mousePosition[0] = 0.5*(
                window_start+window_end)
            self.dataUI.animationPicking.changedSelect = True
            QMessageBox.warning(
                self, 'Manual anchor required',
                f'{len(manual_anchor_requests)} profile branch(es) were not '
                'automatically picked because their closest receiver could not '
                'be validated.\n\n'
                f'The first case is {file_name}, geophone {trace_id+1}. '
                'Pick its first arrival manually (click or A), then run '
                'Automated picking again and choose Yes to preserve it as an anchor.')
        self.statusBar.showMessage(
            f'Hybrid picking finished: {accepted} accepted, {rejected} rejected, '
            f'{defective_count} defective geophone(s) excluded, '
            f'{zero_offset_count} zero-offset t0 reference(s), '
            f'{len(manual_anchor_requests)} branch(es) need a manual anchor. '
            'Rejected traces should be reviewed manually.', 10000)

    def trace_branches(self, file_id):
        """Order receiver traces away from the source on each profile side.

        Parameters
        ----------
        file_id : int
            Index of the shot in ``self.dataUI.sisData``.

        Returns
        -------
        list of numpy.ndarray
            One array of trace identifiers per profile side. Each array starts
            at the receiver closest to the source and progresses outwards.

        Notes
        -----
        For multidimensional coordinates, the profile direction is estimated
        from the receivers by singular-value decomposition.
        """
        nb_traces = len(self.dataUI.sisData[file_id])
        receivers = np.asarray(self.dataUI.geometry.receivers, dtype=float)
        sensors = np.asarray(self.dataUI.geometry.sensors, dtype=float)
        if receivers.ndim == 1:
            receivers = receivers[:, None]
        if sensors.ndim == 1:
            sensors = sensors[:, None]
        if receivers.shape[0] != nb_traces:
            return [np.arange(nb_traces, dtype=int)]

        source_id = int(self.dataUI.geometry.sourcesId[file_id])
        source = sensors[source_id]
        if receivers.shape[1] == 1:
            receiver_position = receivers[:, 0]
            source_position = float(source[0])
        else:
            centred = receivers-np.mean(receivers, axis=0)
            try:
                _, _, axes = np.linalg.svd(centred, full_matrices=False)
                profile_axis = axes[0]
            except np.linalg.LinAlgError:
                profile_axis = np.zeros(receivers.shape[1])
                profile_axis[0] = 1.0
            receiver_position = receivers @ profile_axis
            source_position = float(source @ profile_axis)

        signed_distance = receiver_position-source_position
        left = np.flatnonzero(signed_distance < 0)
        right = np.flatnonzero(signed_distance >= 0)
        branches = []
        for branch in (left, right):
            if branch.size:
                order = branch[np.argsort(np.abs(signed_distance[branch]))]
                branches.append(order.astype(int))
        return branches if branches else [np.arange(nb_traces, dtype=int)]

    def traceNumberChanged(self, value):
        """Select a trace from the spin box and request a plot refresh."""
        self.dataUI.animationPicking.currSelect = value
        curr_pick = self.dataUI.picking[self.dataUI.sisFileId, value]
        self.dataUI.animationPicking.autoPickCenter = (
            None if np.isnan(curr_pick) else float(curr_pick))
        self.dataUI.animationPicking.changedSelect = True

    def updateTab0(self):
        """Enable and refresh picking-tab controls after data-state changes."""
        self.comboBoxFilesPicking.clear()
        for name in self.dataUI.paths.seg2Files:
            self.comboBoxFilesPicking.addItem(name)
        self.dataUI.sisFileId = 0
        self.configure_time_window(reset=True)
        self.comboBoxFilesPicking.currentIndexChanged.connect(
            self.comboBoxChange)
        self.connectMouse = self.mainGraph.mpl_connect(
            'motion_notify_event', self.changeMouse)
        self.connectPress = self.mainGraph.mpl_connect(
            'button_press_event', self.onPress)
        self.connectRelease = self.mainGraph.mpl_connect(
            'button_release_event', self.onRelease)
        self.connectKeyPress = self.mainGraph.mpl_connect(
            'key_press_event', self.onKeyPress)
        if hasattr(self, 'connectScroll'):
            self.mainGraph.mpl_disconnect(self.connectScroll)
        self.connectScroll = self.mainGraph.mpl_connect(
            'scroll_event', self.change_visual_amplitude)
        self.aniMain = animation.FuncAnimation(self.mainGraph.fig, self.animationMain, interval=33.0,
                                               cache_frame_data=False)
        self.aniZoom = animation.FuncAnimation(self.zoomGraph.fig, self.animationZoom, interval=33.0,
                                               cache_frame_data=False)
        self.mainGraph.draw()
        self.zoomGraph.draw()

    def comboBoxChange(self, newId):
        """Switch the active seismic shot and reset its display state."""
        if newId < 0 or newId >= len(self.dataUI.sisData):
            return
        self.dataUI.sisFileId = newId
        self.dataUI.animationPicking.selectedPicks.clear()
        self.configure_time_window(reset=False)
        if self.dataUI.animationPicking.fftShowed:
            self.dataUI.animationPicking.fftShowed = False
            self.fftGraph()
        else:
            self.dataUI.animationPicking.changedSelect = True

    def _savePicking(self):
        """Export valid first-arrival picks and uncertainties to an SGT file."""
        self.statusBar.showMessage('Save current picking . . .')
        # Building the sgt array:
        sensors = self.dataUI.geometry.sensors
        sourcesId = self.dataUI.geometry.sourcesId
        receivers = self.dataUI.geometry.receivers
        picksSave = []
        for nbFile in range(len(self.dataUI.paths.seg2Files)):
            for i in range(len(self.dataUI.sisData[0])):
                sId = int(sourcesId[nbFile])
                rId = int(sensors.index(receivers[i]))
                if sId != rId:  # The traveltime for source = receiever is 0 and not usefull for inversion!
                    t = self.dataUI.picking[nbFile, i]
                    if np.isfinite(t):
                        # Internal errors are relative; pyGIMLi expects an
                        # absolute, strictly positive uncertainty in seconds.
                        dt = float(self.dataUI.sisData[nbFile][0].stats.delta)
                        err = absolute_pick_error(
                            self.dataUI.pickingError[nbFile, i], t, dt)
                        picksSave.append([sId, rId, t, err])
        # Remove unused sensors from the list:
        usedSensors = [False]*len(sensors)
        for pick in picksSave:
            usedSensors[pick[0]] = True
            usedSensors[pick[1]] = True
        oldId = range(len(sensors))
        oldId = [i for i in range(len(sensors)) if usedSensors[i]]
        sensors = [sensors[i] for i in range(len(sensors)) if usedSensors[i]]
        newId = range(len(sensors))
        for pick in picksSave:
            pick[0] = newId[oldId.index(pick[0])]
            pick[1] = newId[oldId.index(pick[1])]
        # Saving the file:
        # The first argument returned is the filename and path
        fname, _ = QFileDialog.getSaveFileName(
            self, 'Select file to save', filter='Source-Receiver-Time file (*.sgt)')
        if fname != "":
            f = open(os.path.join(fname), 'w')  # Create a new file
            nbSensors = len(sensors)
            f.write('%d # shot/geophone points\n' % nbSensors)
            f.write('#x\ty\n')
            for i in range(nbSensors):
                f.write('%.2f\t%.2f\n' % (sensors[i][0], sensors[i][1]))
            nbMeas = len(picksSave)
            f.write('%d # measurements\n' % nbMeas)
            f.write('#s\tg\tt\terr\n')
            for i in range(nbMeas):
                f.write('%d\t%d\t%f\t%f\n' % (
                    picksSave[i][0]+1, picksSave[i][1]+1,
                    picksSave[i][2], picksSave[i][3]))
            f.close()
            self.statusBar.showMessage(f'Picking saved at: {fname}.', 10000)
            self.filePicksPath.setText(fname)
            self._initPygimli(fname)
        else:
            self.statusBar.showMessage('No file saved!', 2000)

    def _loadPicking(self):
        """Load an SGT pick file and initialize the modelling workflow."""
        self.statusBar.showMessage('Loading picking file . . .')
        # 1) Load a file with the first arrival:
        fName, _ = QFileDialog.getOpenFileName(
            self, 'Select file to load', filter='Source-Receiver-Time file (*.sgt)')
        # We retreived a first-arrival file --> geometry of the sensors + first arrivals
        if fName != "":
            with open(fName) as f:
                lines = f.read().splitlines()
            markerNbSensors = "# shot/geophone points"
            markerMeasurements = "# measurements"
            for line in lines:
                if line.endswith(markerNbSensors):
                    nbSensors = int(line[:-len(markerNbSensors)])
                    sensors = np.zeros((nbSensors, 2))
                    idxSensor = 0
                elif line.endswith("#x\ty"):
                    pass
                elif idxSensor < nbSensors:
                    sensors[idxSensor, :] = re.split(r'\t+', line)
                    idxSensor += 1
                elif line.endswith(markerMeasurements):
                    nbMeasurements = int(line[:-len(markerMeasurements)])
                    measurements = np.zeros((nbMeasurements, 4))  # s g t err
                    idxMeas = 0
                elif line.endswith('#s\tg\tt\terr'):
                    measurements = np.zeros((nbMeasurements, 4))  # s g t err
                elif line.endswith('#s\tg\tt'):
                    measurements = np.zeros((nbMeasurements, 3))  # s g t
                elif idxMeas < nbMeasurements:
                    measurements[idxMeas, :] = re.split(r'\t+', line)
                    idxMeas += 1
            self.dataUI.modellingData.sensors = sensors
            self.dataUI.modellingData.measurements = measurements
            loadPicking = False
            if self.dataUI.dataLoaded:
                # Check if the sgt file corresponds to the already loaded dataset (similar array)
                sensorsInList = True
                for i in range(nbSensors):
                    currSensor = sensors[i, :]
                    equals = np.all(
                        self.dataUI.geometry.sensors == currSensor, axis=1)
                    if np.sum(equals) < 1:
                        # The sensor is not in the list
                        sensorsInList = False
                        break
                if sensorsInList:
                    if not (np.all(np.isnan(self.dataUI.picking))):
                        # Load the picking to add atop the traces:
                        reply = QMessageBox.question(
                            self, 'Overwritting picking . . .', 'Are you sure you want to overwrite the current picking?', QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
                        if reply == QMessageBox.Yes:
                            # loading the picking:
                            loadPicking = True
                        else:
                            # Do nothing for the picking window:
                            loadPicking = False
                            self.statusBar.showMessage(
                                f'Data loaded but picking not presented', 2000)
                    else:
                        loadPicking = True
            if loadPicking:
                sources = np.asarray(self.dataUI.geometry.sensors)
                sources = sources[list(
                    self.dataUI.geometry.sourcesId.astype(int)), :]
                receivers = np.asarray(self.dataUI.geometry.receivers)
                self.dataUI.pickingError[:] = np.nan
                self.dataUI.picking[:] = np.nan
                if len(measurements[0, :]) > 3:
                    if np.mean(measurements[:, 3]) < 0.001:
                        reply = QMessageBox.question(
                            self, 'Error scaling issue . . .', 'The *.sgt file you are trying to load has very low mean errors.\nIt is probably an old file with an error in the encoding (<v0.4.0).\nWould you like to overwrite the error with a constant 3%% error ?', QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
                        if reply == QMessageBox.Yes:
                            measurements[:, 3] = 0.03
                for i in range(nbMeasurements):
                    sCurr = sensors[int(measurements[i, 0])-1, :]
                    rCurr = sensors[int(measurements[i, 1])-1, :]
                    pickCurr = measurements[i, 2]
                    if len(measurements[i, :]) > 3:
                        errCurr = measurements[i, 3]
                    sId = np.where(np.all(sources == sCurr, axis=1))[0]
                    rId = np.where(np.all(receivers == rCurr, axis=1))[0]
                    self.dataUI.picking[sId, rId] = pickCurr
                    self.ensure_manual_picking_mask()
                    self.dataUI.manualPickingMask[sId, rId] = True
                    if len(measurements[i, :]) > 3:
                        # /pickCurr # Attention, breaks backward compatibility!!!
                        self.dataUI.pickingError[sId, rId] = errCurr
                    else:
                        self.dataUI.pickingError[sId, rId] = DEFAULT_ERROR
                self.statusBar.showMessage(
                    f'Data loaded with picking on graphs', 10000)
                self.dataUI.animationPicking.changedSelect = True
            self.filePicksPath.setText(fName)
            self._initPygimli(fName)
            self._initModelling()

    def _initPygimli(self, fname):
        """Read an SGT file into pyGIMLi and prepare inversion data objects."""
        # Preparing inversion of data (pygimli)
        self.dataUI.invData.data = pg.DataContainer(fname, sensorTokens='s g')
        self.dataUI.invData.data.sortSensorsX(incX=True)
        self.dataUI.invData.manager = TTMgr(self.dataUI.invData.data)
        self.dataUI.invData.mesh = self._create_inversion_mesh()
        self.invModelGraph.axes.clear()
        self.dataGraph.axes.clear()
        self.fitGraph.axes.clear()
        pgshow(self.dataUI.invData.mesh, ax=self.invModelGraph.axes)
        drawFirstPicks(ax=self.dataGraph.axes, data=self.dataUI.invData.data)
        self.dataGraph.fig.tight_layout()
        self.dataGraph.draw()
        self.invModelGraph.fig.tight_layout()
        self.invModelGraph.draw()

    def _initModelling(self):
        """Build initial layered models for every source and orientation."""
        if len(self.dataUI.modellingData.sensors) == 0:
            return
        sensors, measurements = self.dataUI.modellingData.sensors, self.dataUI.modellingData.measurements
        if np.all(np.abs(sensors[:, -1]) < 1e-6):
            # Plotting the hodochrones:
            self.plotHodochrones(sensors, measurements)
        else:
            QMessageBox.warning(
                self, 'Warning !', 'Impossible to model using the intercept time method!\nNo topography authorized.')
            return
        # Initialize the modelling:
        try:
            self.dataUI.modellingData.nbLayers = self.nbLayersSelector.value()
            self.dataUI.modellingData.hodoPoints = []
            self.dataUI.modellingData.combinationSR = []
            self.dataUI.modellingData.sourcesX = []
            self.dataUI.modellingData.appVelocities = []
            self.dataUI.modellingData.interceptTime = []
            self.dataUI.modellingData.orientations = []
            # Gather the possible sources and orientations:
            self.dataUI.modellingAnimation.namesSources = []
            self.dataUI.modellingAnimation.namesOrientations = []
            sources = np.unique(measurements[:, 0]).astype(int)
            for sId in sources:
                sourceX = sensors[sId-1, 0]
                self.dataUI.modellingAnimation.namesSources.append(
                    f'Source at {sourceX} m.')
                self.dataUI.modellingData.sourcesX.append(sourceX)
                index = measurements[:, 0].astype(int) == sId
                sourceX = sensors[sId-1, 0]
                receiversX = sensors[measurements[index, 1].astype(int) - 1, 0]
                times = measurements[index, 2]
                receiversLeft = receiversX[receiversX < sourceX]
                receiversRight = receiversX[receiversX > sourceX]
                orientationsText = []
                orientations = []
                appVel = []
                intercept = []
                hodoPts = []
                if receiversLeft.size != 0:
                    orientationsText.append('Left')
                    orientations.append(-1)
                    self.dataUI.modellingData.combinationSR.append(
                        [sourceX, -1])
                    try:
                        inter, v, points = build_model(
                            sourceX, receiversLeft, times[receiversX < sourceX], self.dataUI.modellingData.nbLayers, -1)
                        appVel.append(v)
                        intercept.append(inter)
                        hodoPts.append(points)
                    except:
                        QMessageBox.warning(
                            self, 'Warning !', 'Impossible to automatically model!')
                        nbLayers = self.dataUI.modellingData.nbLayers
                        orientation = -1
                        appVel.append(np.linspace(600, 2000, nbLayers))
                        intercept.append(np.linspace(0.0, 0.015, nbLayers))
                        points = np.zeros((nbLayers+1, 2))
                        points[0, 0] = sourceX
                        for i in range(nbLayers):
                            if i < nbLayers-1:
                                xTemp = (inter[i]-inter[i+1]) / \
                                    ((1/v[i+1])-(1/v[i]))
                                tTemp = inter[i] + xTemp*(1/v[i])
                                points[i+1, 0] = sourceX + orientation*xTemp
                                points[i+1, 1] = tTemp
                            else:
                                if orientation > 0:
                                    points[i+1, 0] = max(receiversX)
                                else:
                                    points[i+1, 0] = min(receiversX)
                                points[i+1, 1] = inter[-1] + \
                                    np.abs(sourceX - points[i+1, 0])*(1/v[-1])
                        hodoPts.append(points)
                if receiversRight.size != 0:
                    orientationsText.append('Right')
                    orientations.append(1)
                    self.dataUI.modellingData.combinationSR.append(
                        [sourceX, 1])
                    try:
                        inter, v, points = build_model(
                            sourceX, receiversRight, times[receiversX > sourceX], self.dataUI.modellingData.nbLayers, 1)
                        appVel.append(v)
                        intercept.append(inter)
                        hodoPts.append(points)
                    except:
                        QMessageBox.warning(
                            self, 'Warning !', 'Impossible to automatically model!')
                        nbLayers = self.dataUI.modellingData.nbLayers
                        orientation = 1
                        appVel.append(np.linspace(600, 2000, nbLayers))
                        intercept.append(np.linspace(0.0, 0.015, nbLayers))
                        points = np.zeros((nbLayers+1, 2))
                        points[0, 0] = sourceX
                        for i in range(nbLayers):
                            if i < nbLayers-1:
                                xTemp = (inter[i]-inter[i+1]) / \
                                    ((1/v[i+1])-(1/v[i]))
                                tTemp = inter[i] + xTemp*(1/v[i])
                                points[i+1, 0] = sourceX + orientation*xTemp
                                points[i+1, 1] = tTemp
                            else:
                                if orientation > 0:
                                    points[i+1, 0] = max(receiversX)
                                else:
                                    points[i+1, 0] = min(receiversX)
                                points[i+1, 1] = inter[-1] + \
                                    np.abs(sourceX - points[i+1, 0])*(1/v[-1])
                        hodoPts.append(points)
                self.dataUI.modellingAnimation.namesOrientations.append(
                    orientationsText)
                self.dataUI.modellingData.appVelocities.append(appVel)
                self.dataUI.modellingData.interceptTime.append(intercept)
                self.dataUI.modellingData.orientations.append(orientations)
                self.dataUI.modellingData.hodoPoints.append(hodoPts)
            self._updateHodoGraph()
            self._updateModelGraph()
            # Implement the source/receiver selector:
            self.sourceSelector.clear()
            self.receiversSelector.clear()
            self.sourceSelector.addItems(
                self.dataUI.modellingAnimation.namesSources)
            self.receiversSelector.addItems(
                self.dataUI.modellingAnimation.namesOrientations[0])
            self.sourceSelector.currentIndexChanged.connect(
                self.sourceSelectorChanged)
        except:
            QMessageBox.warning(
                self, 'Warning !', 'Impossible to model using the intercept time method!\nThe robustness of this method needs to be improved.')

    def _updateHodoGraph(self):
        """Redraw measured and interpreted travel-time hodographs."""
        axHod = self.hodochronesGraph.axes
        axHod.cla()
        maxY = self.plotHodochrones(
            self.dataUI.modellingData.sensors, self.dataUI.modellingData.measurements)
        colors = matplotlib.pyplot.cm.tab10(np.arange(10))
        xAxisShow = np.linspace(np.min(self.dataUI.modellingData.sensors[:, 0]), np.max(
            self.dataUI.modellingData.sensors[:, 0]), 1000)
        for i, sourceX in enumerate(self.dataUI.modellingData.sourcesX):
            for j, orientation in enumerate(self.dataUI.modellingData.orientations[i]):
                vel = self.dataUI.modellingData.appVelocities[i][j]
                inter = self.dataUI.modellingData.interceptTime[i][j]
                points = self.dataUI.modellingData.hodoPoints[i][j]
                xShow = xAxisShow[(xAxisShow - sourceX)*orientation >= 0]
                maxY = max([maxY, max(points[:, 1])])
                for k in range(self.dataUI.modellingData.nbLayers):
                    times = inter[k] + (1/vel[k])*np.abs(xShow-sourceX)
                    axHod.plot(xShow, times, color=colors[i % 10])
                axHod.plot(points[:, 0], points[:, 1], linestyle='none',
                           marker='o', color='k', markersize=5)
        axHod.set_ylim(bottom=0.0, top=maxY*1.1)
        self.hodochronesGraph.draw()

    def _updateModelGraph(self):
        """Redraw the layered velocity model derived from the hodographs."""
        axMod = self.modelGraph.axes
        axMod.cla()
        # Show the sensors array at the surface:
        axMod.plot(self.dataUI.modellingData.sensors[:, 0],
                   self.dataUI.modellingData.sensors[:, 1], marker='v', color='k')
        axMod.grid()
        offsetVelocity = self.dataUI.modellingAnimation.offsetVelocity
        sourcesOrientation = []
        for i, sourceX in enumerate(self.dataUI.modellingData.sourcesX):
            for j, orientation in enumerate(self.dataUI.modellingData.orientations[i]):
                sourcesOrientation.append([i, j, sourceX, orientation])
                if orientation > 0:
                    ha = 'left'
                else:
                    ha = 'right'
                vel = self.dataUI.modellingData.appVelocities[i][j]
                inter = self.dataUI.modellingData.interceptTime[i][j]
                _, v, pos = model1D(inter, vel)
                axMod.plot(sourceX + pos[:, 0]*orientation, pos[:, 1],
                           linestyle='none', color='k', marker='x')
                for k in range(self.dataUI.modellingData.nbLayers):
                    axMod.text(sourceX+pos[j, 0] + offsetVelocity, pos[k, 1] +
                               2*offsetVelocity, f'$v_{k}$ = {round(v[k],2)} m/s', ha=ha)
        sourcesOrientation = np.asarray(sourcesOrientation)
        sourcesOrientation[np.lexsort(
            (sourcesOrientation[:, 3], sourcesOrientation[:, 2]))]
        sourcesLeft = sourcesOrientation[sourcesOrientation[:, 3] < 0, :]
        sourcesRight = sourcesOrientation[sourcesOrientation[:, 3] > 0, :]
        for s1 in sourcesRight.tolist():
            sourceX = s1[2]
            for s2 in sourcesLeft[sourcesLeft[:, 2] > sourceX].tolist():
                vel1 = self.dataUI.modellingData.appVelocities[int(
                    s1[0])][int(s1[1])]
                inter1 = self.dataUI.modellingData.interceptTime[int(
                    s1[0])][int(s1[1])]
                vel2 = self.dataUI.modellingData.appVelocities[int(
                    s2[0])][int(s2[1])]
                inter2 = self.dataUI.modellingData.interceptTime[int(
                    s2[0])][int(s2[1])]
                vS = np.asarray([vel1, vel2]).T
                interpS = np.asarray([inter1, inter2]).T
                vSlope, hL, hR = modelWithSlope(interpS, vS)
                for i in range(self.dataUI.modellingData.nbLayers):
                    if i != self.dataUI.modellingData.nbLayers-1:
                        axMod.plot([s1[2], s2[2]], [hL[i], hR[i]], color='k')
                    if i == 0:
                        axMod.text(np.mean([s1[2], s2[2]]), 2*offsetVelocity,
                                   f'$v_{i}$ = {round(vSlope[i],2)} m/s', ha='center')
                    else:
                        axMod.text(np.mean([s1[2], s2[2]]), (hL[i-1] + hR[i-1])/2 + 2 *
                                   offsetVelocity, f'$v_{i}$ = {round(vSlope[i],2)} m/s', ha='center')
        axMod.set_xlabel('X [m]')
        axMod.set_ylabel('Depth [m]')
        xRange = axMod.get_xlim()
        xRange = xRange[1] - xRange[0]
        axMod.set_ylim((-1, np.ceil(xRange/5)))
        axMod.invert_yaxis()
        self.modelGraph.draw()

    def _updateModelling(self, frame):
        """Animation callback that recomputes plots after point edits."""
        # Gather the number of layers in the model
        nbLayers = self.dataUI.modellingData.nbLayers
        # Gathering the info about intercept times and apparent velocities:
        sourceId = self.sourceSelector.currentIndex()
        receiversOrientation = self.receiversSelector.currentIndex()
        orientation = self.dataUI.modellingData.orientations[sourceId][receiversOrientation]
        points = np.asarray(
            self.dataUI.modellingData.hodoPoints[sourceId][receiversOrientation])
        nbLayers = self.dataUI.modellingData.nbLayers
        vel = np.zeros((nbLayers,))
        inter = np.zeros((nbLayers,))
        sourceX = points[0, 0]
        for i in range(self.dataUI.modellingData.nbLayers):
            vel[i] = np.abs((points[i+1, 0] - points[i, 0]) /
                            (points[i+1, 1] - points[i, 1]))
            inter[i] = 1/vel[i] * (sourceX-points[i, 0]) * \
                orientation + points[i, 1]
        self.dataUI.modellingData.appVelocities[sourceId][receiversOrientation] = vel
        self.dataUI.modellingData.interceptTime[sourceId][receiversOrientation] = inter
        # Updating the graphs:
        self._updateHodoGraph()
        self._updateModelGraph()

    def animateModelling(self):
        """Start periodic redraws for interactive layered-model editing."""
        self.connectMouse = self.hodochronesGraph.mpl_connect(
            'motion_notify_event', self.changeMouseModelling)
        self.connectPress = self.hodochronesGraph.mpl_connect(
            'button_press_event', self.onPressModelling)
        self.connectRelease = self.hodochronesGraph.mpl_connect(
            'button_release_event', self.onReleaseModelling)
        self.aniModelling = animation.FuncAnimation(
            self.hodochronesGraph.fig, self._updateModelling)  # , interval=16.7)
        self.modelGraph.draw()
        self.hodochronesGraph.draw()

    def setAnimationModelling(self):
        """Enable or disable mouse-driven hodograph point editing."""
        if self.aniModelling is None:
            if self.movePoints.isChecked():
                self.animateModelling()
        else:
            if self.movePoints.isChecked():
                self.aniModelling.event_source.start()
                self.connectMouse = self.hodochronesGraph.mpl_connect(
                    'motion_notify_event', self.changeMouseModelling)
                self.connectPress = self.hodochronesGraph.mpl_connect(
                    'button_press_event', self.onPressModelling)
                self.connectRelease = self.hodochronesGraph.mpl_connect(
                    'button_release_event', self.onReleaseModelling)
            else:
                self.aniModelling.event_source.stop()
                self.hodochronesGraph.mpl_disconnect(self.connectMouse)
                self.hodochronesGraph.mpl_disconnect(self.connectPress)
                self.hodochronesGraph.mpl_disconnect(self.connectRelease)

    def sourceSelectorChanged(self, i):
        """Update orientation choices and plots for the selected source."""
        self.receiversSelector.clear()
        self.receiversSelector.addItems(
            self.dataUI.modellingAnimation.namesOrientations[i])

    def plotHodochrones(self, sensors, measurements):
        """Plot travel-time measurements grouped by source location."""
        ax = self.hodochronesGraph.axes
        maxY = 0
        sources = np.unique(measurements[:, 0]).astype(int)
        colors = matplotlib.pyplot.cm.tab10(np.arange(10))
        hasError = len(measurements[0, :]) == 4
        for i, sId in enumerate(sources):
            index = measurements[:, 0].astype(int) == sId
            ti = measurements[index, 2]
            maxY = max([maxY, max(ti)])
            ri = sensors[measurements[index, 1].astype(int) - 1, 0]
            if hasError:
                erri = measurements[index, 3]
                ax.errorbar(ri, ti, yerr=np.multiply(
                    erri, ti), linestyle='none', color=colors[i % 10], marker='s', markersize=5)
            else:
                ax.plot(ri, ti, linestyle='none',
                        color=colors[i % 10], marker='s', markersize=5)
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Traveltime (s)')
        ax.set_ylim(bottom=0.0)
        ax.grid(True)
        self.hodochronesGraph.fig.tight_layout()
        self.hodochronesGraph.draw()
        return maxY

    # def _saveModel(self):
    #     # TODO
    #     self.statusBar.showMessage('Saving current model . . .')
    #     time.sleep(10)
    #     self.statusBar.showMessage(DEFAULT_STATUS)

    def _loadInvMesh(self):
        """Load a previously saved mesh for the inversion."""
        # Load an inversion mesh that was already created (*.poly)
        fName, _ = QFileDialog.getOpenFileName(
            self, 'Select file to load', filter='GIMLi mesh file (*.bms)')
        if fName != "":
            self.dataUI.invData.mesh = pg.Mesh(fName)
            self.invModelGraph.axes.cla()
            pg.show(self.dataUI.invData.mesh, ax=self.invModelGraph.axes)
            self.invModelGraph.draw()
            self.dataUI.meshLoaded = True
            self.statusBar.showMessage(f'Mesh loaded from {fName}', 10000)
        else:
            self.statusBar.showMessage('No mesh could be loaded!', 2000)

    def _loadInitModel(self):
        """Load a starting velocity model for the active inversion mesh."""
        # Load an existing model as the starting model (*.vector) for the inversion
        fName, _ = QFileDialog.getOpenFileName(
            self, 'Select file to load', filter='Result Vector file (*.vector)')
        if fName != "":
            self.dataUI.invData.startModel = pg.Vector(np.loadtxt(fName))
            self.invModelGraph.axes.cla()
            try:
                pg.show(self.dataUI.invData.mesh,
                        self.dataUI.invData.startModel, ax=self.invModelGraph.axes)
                self.invModelGraph.draw()
                self.statusBar.showMessage(
                    f'Initial model loaded from {fName}', 10000)
            except:
                QMessageBox.warning(
                    self, 'Warning !', 'The initial model and the current mesh do not match in size!')
                self.dataUI.invData.startModel = None
                self.statusBar.showMessage('No initial model loaded!', 2000)
        else:
            self.statusBar.showMessage('No initial model loaded!', 2000)

    def _saveInvMesh(self):
        """Save the current inversion mesh to a user-selected file."""
        # Save the inversion mesh (*.poly)
        fName, _ = QFileDialog.getSaveFileName(
            self, 'Select file to save', filter='GIMLi mesh file (*.bms)')
        if fName != "":
            self.dataUI.invData.mesh.save(fName)
            self.statusBar.showMessage(f'Mesh saved to {fName}', 10000)
        else:
            self.statusBar.showMessage('No mesh saved!', 2000)

    def _saveInvAsVTK(self):
        """Export the inversion mesh and model values in VTK format."""
        # Save the inversion results into a VTK file (for Paraview)
        fName, _ = QFileDialog.getSaveFileName(
            self, 'Select file to save', filter='Paraview mesh file (*.vtk)')
        if fName != "":
            mgr = self.dataUI.invData.manager
            m = self.dataUI.invData.mesh
            m.addData("Velocity [m/s]", mgr.model)
            coverage = mgr.fop.jacobian().transMult(np.ones(mgr.fop.data.size()))
            m.addData("Coverage [/]", coverage)
            C = mgr.fop.constraintsRef()
            m.addData(
                "Standardized Coverage [/]", np.sign(np.absolute(C.transMult(C * coverage))))
            m.exportVTK(fName)
            self.statusBar.showMessage(f'Results saved to {fName}', 10000)
        else:
            self.statusBar.showMessage('Results where NOT saved!', 2000)

    def _saveInvResponse(self):
        """Export the simulated travel-time response of the inversion result."""
        # Save the model response for the last iteration (*.vector)
        fName, _ = QFileDialog.getSaveFileName(
            self, 'Select file to save', filter='Response Vector file (*.vector)')
        if fName != "":
            np.savetxt(fName, self.dataUI.invData.manager.inv.response)
            self.statusBar.showMessage(f'Response saved to {fName}', 10000)
        else:
            self.statusBar.showMessage('Response was NOT saved!', 2000)

    def _saveInvResult(self):
        """Save the recovered velocity model as a numeric vector."""
        # Save the model for the last iteration (*.vector)
        fName, _ = QFileDialog.getSaveFileName(
            self, 'Select file to save', filter='Result Vector file (*.vector)')
        if fName != "":
            np.savetxt(fName, self.dataUI.invData.manager.model)
            self.statusBar.showMessage(f'Result saved to {fName}', 10000)
        else:
            self.statusBar.showMessage('Result was NOT saved!', 2000)

    def _pickTracesTabUI(self):
        """Construct and return the seismic trace-picking tab."""
        importTab = QWidget(self.tabs)
        layout = QGridLayout(self.tabs)  # Grid of 10-by-15

        # Comb Box for choosing the correct file to pick.
        self.comboBoxFilesPicking = QComboBox()
        for name in self.dataUI.paths.seg2Files:
            self.comboBoxFilesPicking.addItem(name)
        layout.addWidget(self.comboBoxFilesPicking, 0, 0, 1, 15)

        # Main graph with all the traces
        self.mainGraph = MplCanvas(importTab, width=8, height=7, dpi=100)
        self.mainGraph.axes.plot(
            [0, 1, 2, 3, 4], [10, 1, 20, 3, 40], animated=True)
        self.mainGraphToolbar = CustomHomeToolbar(self.mainGraph, importTab)
        mainGraphLayout = QVBoxLayout()
        mainGraphLayout.addWidget(self.mainGraphToolbar)
        mainGraphLayout.addWidget(self.mainGraph)
        layout.addLayout(mainGraphLayout, 2, 0, 9, 10)
        # Zoom graph with the current trace
        self.zoomGraph = MplCanvas(importTab, width=5, height=4, dpi=100)
        self.zoomGraph.axes.plot(
            [0, 1, 2, 3, 4], [10, 1, 20, 3, 40], animated=True)
        self.zoomGraph.axes.tick_params(
            axis='both',       # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            left=False,
            right=False,
            labelleft=False,
            labelbottom=False)  # labels along the bottom edge are off
        self.aniMain = None
        self.aniZoom = None
        layout.addWidget(self.zoomGraph, 2, 10, 4, 5)
        self.buttonTabPickingSet = QPushButton('Set picking')
        self.buttonTabPickingReset = QPushButton('Reset picking')
        self.buttonManageDefective = QPushButton('Trace validity…')
        self.buttonManageDefective.setToolTip(
            'Mark the selected trace as defective or restore it as valid, '
            'for this shot or all shots.')
        self.buttonManageDefective.setEnabled(False)
        self.spinBoxCurrSelect = QSpinBox(importTab)
        self.spinBoxCurrSelect.setEnabled(False)
        self.buttonTabPickingSetT0 = QPushButton('Set t=0')
        self.buttonTabPickingSetT0.clicked.connect(self.setT0)
        textBox = QLabel(importTab)
        textBox.setText('Select the trace number:')
        layout.addWidget(textBox, 6, 10, 1, 5)
        layout.addWidget(self.spinBoxCurrSelect, 7, 10, 1, 5)
        layout.addWidget(self.buttonTabPickingSetT0, 8, 10, 1, 5)
        layout.addWidget(self.buttonTabPickingSet, 9, 10, 1, 5)
        layout.addWidget(self.buttonTabPickingReset, 10, 10, 1, 5)
        layout.addWidget(self.buttonManageDefective, 11, 10, 1, 5)
        shortcutHelp = QLabel(
            'Shortcuts: Up = trace above, Down = trace below  |  Left/Right fine adjust\n'
            'Shift+Left/Right: 10 samples  |  A: automatic pick\n'
            'Drag: select picks  |  Ctrl+Drag: line picking\n'
            'Delete/Backspace: remove selected picks (blue)\n'
            'Automatic picks: green = reliable, orange = review\n'
            'Grey dashed = defective (automatic or manual); blue dotted = zero-offset t0 reference')
        shortcutHelp.setWordWrap(True)
        shortcutHelp.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop)
        # Keep the right column dedicated to the active trace and controls.
        layout.addWidget(shortcutHelp, 11, 0, 1, 10)
        self.visualAmplitudeLabel = QLabel(
            'Visual amplitude: ×1.00 (mouse wheel over the gather)')
        layout.addWidget(self.visualAmplitudeLabel, 12, 0, 1, 10)

        timeWindowLabel = QLabel('Global signal time window [s]:')
        layout.addWidget(timeWindowLabel, 12, 10, 1, 5)
        self.timeWindowStart = QDoubleSpinBox(importTab)
        self.timeWindowEnd = QDoubleSpinBox(importTab)
        for spinbox in (self.timeWindowStart, self.timeWindowEnd):
            spinbox.setDecimals(6)
            spinbox.setKeyboardTracking(False)
            spinbox.setEnabled(False)
        self.timeWindowStart.setPrefix('From ')
        self.timeWindowEnd.setPrefix('To ')
        timeWindowInputs = QHBoxLayout()
        timeWindowInputs.addWidget(self.timeWindowStart)
        timeWindowInputs.addWidget(self.timeWindowEnd)
        layout.addLayout(timeWindowInputs, 13, 10, 1, 5)

        self.buttonApplyTimeWindow = QPushButton('Apply time window')
        self.buttonResetTimeWindow = QPushButton('Full signal')
        self.buttonApplyTimeWindow.setEnabled(False)
        self.buttonResetTimeWindow.setEnabled(False)
        self.buttonApplyTimeWindow.clicked.connect(self.apply_time_window)
        self.buttonResetTimeWindow.clicked.connect(self.reset_time_window)
        timeWindowButtons = QHBoxLayout()
        timeWindowButtons.addWidget(self.buttonApplyTimeWindow)
        timeWindowButtons.addWidget(self.buttonResetTimeWindow)
        layout.addLayout(timeWindowButtons, 14, 10, 1, 5)

        self.buttonTabPickingSet.setCheckable(True)
        self.buttonTabPickingSet.clicked.connect(self.setPicking)
        self.buttonTabPickingReset.clicked.connect(self.resetPicking)
        self.buttonManageDefective.clicked.connect(
            self.manage_current_defective_trace)
        # self.buttonTabPickingSetT0.setCheckable(True)
        importTab.setLayout(layout)
        return importTab

    def select_t0_near_branch(self, ordered_ids, distances, picks, dt):
        """Keep the closest coherent travel-time branch for t0 estimation.

        Parameters
        ----------
        ordered_ids : array-like of int
            Candidate trace identifiers ordered by increasing source distance.
        distances : array-like
            Source-receiver distances in metres, indexed by trace identifier.
        picks : array-like
            Travel-time picks in seconds, indexed by trace identifier.
        dt : float
            Sampling interval in seconds.

        Returns
        -------
        numpy.ndarray of int
            Identifiers belonging to the initial coherent branch. Callers
            provide at most five candidates.
        """
        ordered_ids = np.asarray(ordered_ids, dtype=int)
        if ordered_ids.size <= 2:
            return ordered_ids

        selected = [int(ordered_ids[0]), int(ordered_ids[1])]
        for trace_id in ordered_ids[2:]:
            x = np.asarray([distances[i] for i in selected], dtype=float)
            y = np.asarray([picks[i] for i in selected], dtype=float)
            if np.ptp(x) <= np.finfo(float).eps:
                continue
            if len(selected) >= 3:
                regression = stats.theilslopes(y, x)
                slope = float(regression.slope if hasattr(regression, 'slope')
                              else regression[0])
                intercept = float(
                    regression.intercept if hasattr(regression, 'intercept')
                    else regression[1])
            else:
                slope, intercept = np.polyfit(x, y, 1)
            predicted = slope*distances[trace_id] + intercept
            residual = abs(picks[trace_id]-predicted)
            existing_residuals = np.abs(y-(slope*x+intercept))
            scatter = (1.4826*np.median(existing_residuals)
                       if existing_residuals.size else 0.0)
            time_tolerance = max(8.0*abs(dt), 7.5e-4, 3.0*scatter)

            previous_id = selected[-1]
            distance_step = distances[trace_id]-distances[previous_id]
            if distance_step <= np.finfo(float).eps:
                continue
            local_slope = (picks[trace_id]-picks[previous_id])/distance_step
            slope_tolerance = max(
                0.50*abs(slope), time_tolerance/distance_step)
            slope_break = abs(local_slope-slope) > slope_tolerance
            if residual > time_tolerance and slope_break:
                break
            selected.append(int(trace_id))
        return np.asarray(selected, dtype=int)

    def recalculateT0FromGeometry(self):
        """Shift every shot so its inferred zero-distance travel time is zero.

        Returns
        -------
        None

        Notes
        -----
        A colocated source-receiver pick is preferred when available. Offset
        shots instead use a robust extrapolation of up to five nearby picks on
        the first coherent branch. After user confirmation, the correction is
        subtracted from finite picks and from ``self.dataUI.beginTime``.
        """
        if not self.dataUI.dataLoaded or len(self.dataUI.sisData) == 0:
            QMessageBox.warning(self, 'Warning!', 'No seismic data are loaded.')
            return

        receivers = np.asarray(self.dataUI.geometry.receivers, dtype=float)
        sensors = np.asarray(self.dataUI.geometry.sensors, dtype=float)
        corrections = np.full(len(self.dataUI.sisData), np.nan)
        methods = []
        for file_id in range(len(self.dataUI.sisData)):
            source_id = int(self.dataUI.geometry.sourcesId[file_id])
            distances = calculateDistance(receivers, sensors[source_id])
            picks = np.asarray(self.dataUI.picking[file_id], dtype=float)
            valid = np.isfinite(picks)
            valid &= ~self.effective_defective_mask(file_id)
            if not np.any(valid):
                methods.append('unavailable')
                continue

            positive_distances = distances[distances > 0]
            zero_tolerance = (max(1e-9, 1e-6*np.min(positive_distances))
                              if positive_distances.size else 1e-9)
            colocated = valid & (distances <= zero_tolerance)
            if np.any(colocated):
                corrections[file_id] = float(np.median(picks[colocated]))
                methods.append('zero-offset receiver')
                continue

            valid_ids = np.flatnonzero(valid)
            valid_ids = valid_ids[np.argsort(distances[valid_ids])]
            # Use at most five nearby receivers, but stop as soon as the
            # closest travel-time branch exhibits a significant slope break.
            # This prevents a deeper/faster phase from biasing the t0
            # extrapolation of an offset shot.
            trial_ids = valid_ids[:min(5, valid_ids.size)]
            dt = float(self.dataUI.sisData[file_id][0].stats.delta)
            near_ids = self.select_t0_near_branch(
                trial_ids, distances, picks, dt)
            near_distances = distances[near_ids]
            near_picks = picks[near_ids]
            if len(np.unique(near_distances)) < 2:
                methods.append('unavailable')
                continue
            if near_ids.size >= 3:
                regression = stats.theilslopes(near_picks, near_distances)
                corrections[file_id] = float(
                    regression.intercept
                    if hasattr(regression, 'intercept') else regression[1])
                methods.append('robust offset extrapolation')
            else:
                _, intercept = np.polyfit(near_distances, near_picks, 1)
                corrections[file_id] = float(intercept)
                methods.append('two-point offset extrapolation')

        available = np.isfinite(corrections)
        if not np.any(available):
            QMessageBox.warning(
                self, 't0 correction',
                'No t0 correction can be estimated from the current picks.')
            return
        reply = QMessageBox.question(
            self, 'Recalculate t=0 from geometry',
            f'A t0 correction can be estimated for {int(np.sum(available))} '
            f'of {len(corrections)} shots.\n'
            f'Correction range: {np.nanmin(corrections):.6f} to '
            f'{np.nanmax(corrections):.6f} s.\n\nApply these corrections?',
            QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
        if reply != QMessageBox.Yes:
            return

        for file_id, correction in enumerate(corrections):
            if not np.isfinite(correction):
                continue
            finite_picks = np.isfinite(self.dataUI.picking[file_id])
            self.dataUI.picking[file_id, finite_picks] -= correction
            self.dataUI.beginTime[file_id] -= correction
        current_correction = corrections[self.dataUI.sisFileId]
        if (np.isfinite(current_correction) and
                self.dataUI.animationPicking.autoPickCenter is not None):
            self.dataUI.animationPicking.autoPickCenter -= current_correction
            self.dataUI.animationPicking.mousePosition[0] -= current_correction
        self.configure_time_window(reset=True)
        self.dataUI.animationPicking.changedSelect = True
        direct_count = sum(method == 'zero-offset receiver' for method in methods)
        self.statusBar.showMessage(
            f't0 corrected for {int(np.sum(available))} shots: '
            f'{direct_count} direct, {int(np.sum(available))-direct_count} extrapolated.',
            10000)

    def setT0(self):
        """Open the t0 helper and apply accepted per-shot time corrections."""
        if self.dataUI.dataLoaded:
            newWindow = PickT0(self)
            newWindow.show()
            newWindow.exec()
            newT0 = newWindow.getNewT0()
            newWindow.close()
            for i, t0Update in enumerate(newT0):
                self.dataUI.beginTime[i] -= t0Update
                self.dataUI.picking[i, :] -= t0Update
            self.dataUI.animationPicking.changedSelect = True

    def setPicking(self):
        """Suspend or resume plot interaction while manual-picking mode changes."""
        if self.dataUI.dataLoaded:
            if self.buttonTabPickingSet.isChecked():
                # We stop the animation:
                self.aniZoom.event_source.stop()
                self.aniMain.event_source.stop()
                self.mainGraph.mpl_disconnect(self.connectMouse)
                self.mainGraph.mpl_disconnect(self.connectPress)
                self.mainGraph.mpl_disconnect(self.connectRelease)
                self.mainGraph.mpl_disconnect(self.connectKeyPress)
                self.spinBoxCurrSelect.setEnabled(False)
                self.buttonManageDefective.setEnabled(False)
            else:
                # We begin back the animation:
                self.aniZoom.event_source.start()
                self.aniMain.event_source.start()
                self.connectMouse = self.mainGraph.mpl_connect(
                    'motion_notify_event', self.changeMouse)
                self.connectPress = self.mainGraph.mpl_connect(
                    'button_press_event', self.onPress)
                self.connectRelease = self.mainGraph.mpl_connect(
                    'button_release_event', self.onRelease)
                self.connectKeyPress = self.mainGraph.mpl_connect(
                    'key_press_event', self.onKeyPress)
                self.spinBoxCurrSelect.setEnabled(True)
                self.buttonManageDefective.setEnabled(True)

    def resetPicking(self):
        """Clear all picks and manual-pick flags from the loaded dataset."""
        if self.dataUI.dataLoaded:
            self.dataUI.picking = np.empty(
                (len(self.dataUI.paths.seg2Files), len(self.dataUI.sisData[0])))
            self.dataUI.picking[:] = np.nan
            self.dataUI.manualPickingMask = np.zeros_like(
                self.dataUI.picking, dtype=bool)
            self.dataUI.animationPicking.selectedPicks.clear()
            self.dataUI.animationPicking.changedSelect = True

    def nearest_hodograph_pick(self, x_value, time_value):
        """Return the shot/trace whose displayed pick is nearest a graph click."""
        picks = np.asarray(self.dataUI.picking, dtype=float)
        receivers = np.asarray(self.dataUI.geometry.receivers, dtype=float)
        if (picks.ndim != 2 or receivers.ndim != 2 or
                receivers.shape[0] != picks.shape[1] or
                not np.isfinite(x_value) or not np.isfinite(time_value)):
            return None

        x_limits = self.dataGraph.axes.get_xlim()
        y_limits = self.dataGraph.axes.get_ylim()
        x_scale = max(abs(x_limits[1]-x_limits[0]), np.finfo(float).eps)
        y_scale = max(abs(y_limits[1]-y_limits[0]), np.finfo(float).eps)
        best = None
        sensors = self.dataUI.geometry.sensors
        source_ids = self.dataUI.geometry.sourcesId
        for file_id in range(min(picks.shape[0], len(source_ids))):
            source_id = int(source_ids[file_id])
            if source_id < 0 or source_id >= len(sensors):
                continue
            source = np.asarray(sensors[source_id], dtype=float)
            valid = np.isfinite(picks[file_id])
            # Source=receiver rows are intentionally absent from the .sgt and
            # from pyGIMLi's hodograph, so do not let them be selected here.
            valid &= np.any(np.abs(receivers-source) > 1e-9, axis=1)
            for trace_id in np.flatnonzero(valid):
                distance = np.hypot(
                    (receivers[trace_id, 0]-x_value)/x_scale,
                    (picks[file_id, trace_id]-time_value)/y_scale)
                if best is None or distance < best[0]:
                    best = (float(distance), file_id, int(trace_id))
        return best

    def nearest_hodograph_sgt_pick(self, x_value, time_value):
        """Identify a hodograph directly from a loaded .sgt file."""
        sensors = np.asarray(self.dataUI.modellingData.sensors, dtype=float)
        measurements = np.asarray(
            self.dataUI.modellingData.measurements, dtype=float)
        if (sensors.ndim != 2 or measurements.ndim != 2 or
                measurements.shape[1] < 3 or not np.isfinite(x_value) or
                not np.isfinite(time_value)):
            return None
        x_limits = self.dataGraph.axes.get_xlim()
        y_limits = self.dataGraph.axes.get_ylim()
        x_scale = max(abs(x_limits[1]-x_limits[0]), np.finfo(float).eps)
        y_scale = max(abs(y_limits[1]-y_limits[0]), np.finfo(float).eps)
        best = None
        for row in measurements:
            source_id, receiver_id = int(row[0]), int(row[1])
            if (receiver_id < 1 or receiver_id > len(sensors) or
                    not np.isfinite(row[2])):
                continue
            distance = np.hypot(
                (sensors[receiver_id-1, 0]-x_value)/x_scale,
                (row[2]-time_value)/y_scale)
            if best is None or distance < best[0]:
                best = (float(distance), source_id)
        return best

    def onHodographClick(self, event):
        """Display the SEG-2 shot name associated with a clicked hodograph."""
        if (event.button != MouseButton.LEFT or
                event.inaxes is not self.dataGraph.axes):
            return
        nearest = self.nearest_hodograph_pick(event.xdata, event.ydata)
        if nearest is not None:
            distance, file_id, trace_id = nearest
            if distance <= 0.06:
                shot_name = os.path.basename(
                    self.dataUI.paths.seg2Files[file_id])
                message = f'Shot: {shot_name}  (trace {trace_id})'
                self.hodographShotLabel.setText(message)
                self.statusBar.showMessage(message, 5000)
                return
        sgt_nearest = self.nearest_hodograph_sgt_pick(
            event.xdata, event.ydata)
        if sgt_nearest is None:
            return
        distance, source_id = sgt_nearest
        # Avoid assigning a shot when the click was plainly in empty space.
        if distance > 0.06:
            self.hodographShotLabel.setText(
                'No hodograph close to this click.')
            return
        sensors = np.asarray(self.dataUI.modellingData.sensors, dtype=float)
        if 1 <= source_id <= len(sensors):
            source = sensors[source_id-1]
            message = (f'Shot sensor {source_id}: '
                       f'({source[0]:.2f}, {source[1]:.2f})')
        else:
            message = f'Shot sensor {source_id}'
        self.hodographShotLabel.setText(message)
        self.statusBar.showMessage(message, 5000)

    def reciprocity_mismatches(self):
        """Find reciprocal picks whose difference exceeds three uncertainties.

        Returns
        -------
        mismatches : list of dict
            Violating shot pairs and their two times, absolute difference,
            three-sigma threshold, and combined one-sigma uncertainty.
        checked : int
            Number of usable reciprocal pairs tested.

        Notes
        -----
        The calculation uses picks currently held by the interface. Offset
        shots, missing picks, and traces marked defective are ignored.
        """
        picks = np.asarray(self.dataUI.picking, dtype=float)
        errors = np.asarray(self.dataUI.pickingError, dtype=float)
        receivers = np.asarray(self.dataUI.geometry.receivers, dtype=float)
        sensors = np.asarray(self.dataUI.geometry.sensors, dtype=float)
        source_ids = np.asarray(self.dataUI.geometry.sourcesId, dtype=int)
        if (picks.ndim != 2 or errors.shape != picks.shape or
                receivers.ndim != 2 or receivers.shape[0] != picks.shape[1]):
            return [], 0

        shot_count = min(picks.shape[0], len(source_ids),
                         len(self.dataUI.sisData))
        source_receivers = []
        for file_id in range(shot_count):
            source_id = source_ids[file_id]
            if source_id < 0 or source_id >= len(sensors):
                source_receivers.append(None)
                continue
            matching = np.flatnonzero(np.all(np.isclose(
                receivers, sensors[source_id], rtol=0.0, atol=1e-7), axis=1))
            source_receivers.append(int(matching[0]) if matching.size else None)

        mismatches = []
        checked = 0
        for first_shot in range(shot_count):
            first_receiver = source_receivers[first_shot]
            if first_receiver is None:
                continue  # Offset shot: no receiver at this source location.
            excluded_first = self.effective_defective_mask(first_shot)
            for second_shot in range(first_shot + 1, shot_count):
                second_receiver = source_receivers[second_shot]
                if second_receiver is None:
                    continue
                excluded_second = self.effective_defective_mask(second_shot)
                if (excluded_first[second_receiver] or
                        excluded_second[first_receiver]):
                    continue
                first_time = picks[first_shot, second_receiver]
                second_time = picks[second_shot, first_receiver]
                if not (np.isfinite(first_time) and np.isfinite(second_time)):
                    continue
                checked += 1
                first_dt = float(
                    self.dataUI.sisData[first_shot][0].stats.delta)
                second_dt = float(
                    self.dataUI.sisData[second_shot][0].stats.delta)
                first_error = absolute_pick_error(
                    errors[first_shot, second_receiver], first_time, first_dt)
                second_error = absolute_pick_error(
                    errors[second_shot, first_receiver], second_time, second_dt)
                combined_error = np.hypot(first_error, second_error)
                difference = abs(first_time-second_time)
                threshold = 3.0*combined_error
                if difference > threshold:
                    mismatches.append({
                        'first_shot': first_shot,
                        'second_shot': second_shot,
                        'first_time': float(first_time),
                        'second_time': float(second_time),
                        'difference': float(difference),
                        'threshold': float(threshold),
                        'sigma': float(combined_error)})
        mismatches.sort(key=lambda item: item['difference']/item['threshold'],
                        reverse=True)
        return mismatches, checked

    def reciprocity_mismatches_from_sgt(self):
        """Check reciprocal pairs directly from the currently loaded SGT data.

        Returns
        -------
        tuple or None
            ``(mismatches, checked)`` with the same meaning as
            :meth:`reciprocity_mismatches`, or ``None`` when the loaded
            modelling arrays do not contain valid SGT measurements.
        """
        sensors = np.asarray(self.dataUI.modellingData.sensors, dtype=float)
        measurements = np.asarray(
            self.dataUI.modellingData.measurements, dtype=float)
        if (sensors.ndim != 2 or measurements.ndim != 2 or
                measurements.shape[1] < 3):
            return None

        records = {}
        for row in measurements:
            source_id, receiver_id = int(row[0]), int(row[1])
            travel_time = float(row[2])
            if source_id == receiver_id or not np.isfinite(travel_time):
                continue
            if measurements.shape[1] >= 4 and np.isfinite(row[3]):
                absolute_error = max(abs(float(row[3])), 1e-6)
            else:
                absolute_error = absolute_pick_error(
                    DEFAULT_ERROR, travel_time, 0.0)
            # A duplicated measurement should not normally occur; keep the
            # most precise one if it does.
            key = (source_id, receiver_id)
            if key not in records or absolute_error < records[key][1]:
                records[key] = (travel_time, absolute_error)

        mismatches = []
        checked = 0
        for (source_id, receiver_id), (first_time, first_error) in records.items():
            reciprocal_key = (receiver_id, source_id)
            if source_id >= receiver_id or reciprocal_key not in records:
                continue
            second_time, second_error = records[reciprocal_key]
            checked += 1
            combined_error = np.hypot(first_error, second_error)
            difference = abs(first_time-second_time)
            threshold = 3.0*combined_error
            if difference > threshold:
                def label(sensor_id):
                    """Format a sensor identifier and its coordinates for reports."""
                    if 1 <= sensor_id <= len(sensors):
                        xy = sensors[sensor_id-1]
                        return f'shot sensor {sensor_id} ({xy[0]:.2f}, {xy[1]:.2f})'
                    return f'shot sensor {sensor_id}'
                mismatches.append({
                    'first_label': label(source_id),
                    'second_label': label(receiver_id),
                    'first_time': first_time,
                    'second_time': second_time,
                    'difference': difference,
                    'threshold': threshold,
                    'sigma': combined_error})
        mismatches.sort(key=lambda item: item['difference']/item['threshold'],
                        reverse=True)
        return mismatches, checked

    def checkReciprocity(self):
        """Report significant violations of seismic reciprocity."""
        from_sgt = self.reciprocity_mismatches_from_sgt()
        if from_sgt is not None:
            mismatches, checked = from_sgt
        elif self.dataUI.dataLoaded and len(self.dataUI.sisData) > 0:
            mismatches, checked = self.reciprocity_mismatches()
        else:
            QMessageBox.warning(self, 'Reciprocity',
                                'Load an .sgt file or seismic data with picks first.')
            return
        if checked == 0:
            QMessageBox.information(
                self, 'Reciprocity',
                'No usable reciprocal shot/receiver pair was found.\n\n'
                'Both shot locations must coincide with receiver locations, '
                'and both picks must be valid.')
            return
        if not mismatches:
            QMessageBox.information(
                self, 'Reciprocity',
                f'{checked} reciprocal pair(s) checked: no significant mismatch '\
                'above 3σ.')
            return

        names = self.dataUI.paths.seg2Files
        lines = []
        for item in mismatches:
            if 'first_label' in item:
                first_name = item['first_label']
                second_name = item['second_label']
            else:
                first_name = os.path.basename(names[item['first_shot']])
                second_name = os.path.basename(names[item['second_shot']])
            lines.append(
                f'{first_name} ↔ {second_name}: '
                f'Δt = {1000*item["difference"]:.2f} ms '
                f'(3σ = {1000*item["threshold"]:.2f} ms; '
                f'{item["first_time"]:.6f} s vs '
                f'{item["second_time"]:.6f} s)')
        message = QMessageBox(self)
        message.setIcon(QMessageBox.Warning)
        message.setWindowTitle('Reciprocity check')
        message.setText(
            f'{len(mismatches)} significant mismatch(es) among {checked} '\
            'reciprocal pair(s).')
        message.setInformativeText(
            'A mismatch is reported when |tAB − tBA| exceeds three combined '\
            'standard uncertainties. This assumes a consistent t0 correction.')
        message.setDetailedText('\n'.join(lines))
        message.exec()

    def _inversionTabUI(self):
        """Construct and return the travel-time inversion configuration tab."""
        # Tab for the inversion of the data using pygimli api.
        inversionTab = QWidget(self.tabs)
        layout = QGridLayout(self.tabs)
        # Selecting the picking file:
        self.filePicksPath = QLabel('File path', inversionTab)
        self.filePicksPath.setAlignment(QtCore.Qt.AlignCenter)
        self.buttonGetSgtFile = QPushButton('...', inversionTab)
        self.buttonGetSgtFile.clicked.connect(self._loadPicking)
        layout.addWidget(self.filePicksPath, 0, 0, 1, 9)
        layout.addWidget(self.buttonGetSgtFile, 0, 9, 1, 1)
        self.invModelGraph = MplCanvas(inversionTab)
        invModelGraphToolbar = NavigationToolbar2QT(
            self.invModelGraph, inversionTab)
        invModelGraphLayout = QVBoxLayout()
        invModelGraphLayout.addWidget(invModelGraphToolbar)
        invModelGraphLayout.addWidget(self.invModelGraph)
        layout.addLayout(invModelGraphLayout, 1, 0, 5, 5)
        self.dataGraph = MplCanvas(inversionTab)
        dataGraphToolbar = NavigationToolbar2QT(self.dataGraph, inversionTab)
        dataGraphLayout = QVBoxLayout()
        dataGraphLayout.addWidget(
            dataGraphToolbar, alignment=QtCore.Qt.AlignRight)
        dataGraphLayout.addWidget(self.dataGraph)
        self.hodographShotLabel = QLabel(
            'Click a hodograph to identify its shot.', inversionTab)
        self.hodographShotLabel.setAlignment(QtCore.Qt.AlignCenter)
        dataGraphLayout.addWidget(self.hodographShotLabel)
        self.buttonCheckReciprocity = QPushButton(
            'Check reciprocity', inversionTab)
        self.buttonCheckReciprocity.setToolTip(
            'Compare reciprocal first arrivals and list significant mismatches.')
        self.buttonCheckReciprocity.clicked.connect(self.checkReciprocity)
        dataGraphLayout.addWidget(self.buttonCheckReciprocity)
        self.connectHodographPick = self.dataGraph.mpl_connect(
            'button_press_event', self.onHodographClick)
        layout.addLayout(dataGraphLayout, 1, 5, 5, 5)
        self.fitGraph = MplCanvas(inversionTab)
        fitGraphToolbar = NavigationToolbar2QT(self.fitGraph, inversionTab)
        fitGraphLayout = QVBoxLayout()
        fitGraphLayout.addWidget(self.fitGraph)
        fitGraphLayout.addWidget(
            fitGraphToolbar, alignment=QtCore.Qt.AlignRight)
        layout.addLayout(fitGraphLayout, 6, 5, 5, 5)
        self.groupeOption = QGroupBox(inversionTab)
        self.groupeOption.setTitle('Inversion options')
        # List of inversion options to enable:
        groupOptionLayout = QGridLayout(self.groupeOption)  # Grid of 5 by 4
        self.setLambda = QLineEdit(
            str(self.dataUI.invData.lam), self.groupeOption)
        self.setLambda.setValidator(
            QtGui.QDoubleValidator(0.0, 10000.0, 2, self.setLambda))
        lambdaText = QLabel('Lambda :')
        self.setZWeight = QLineEdit(
            str(self.dataUI.invData.zWeight), self.groupeOption)
        self.setZWeight.setValidator(
            QtGui.QDoubleValidator(0.01, 100.0, 3, self.setZWeight))
        zWeightText = QLabel('Z-weight (/) :')
        self.setVTop = QLineEdit(
            str(self.dataUI.invData.vTop), self.groupeOption)
        self.setVTop.setValidator(
            QtGui.QDoubleValidator(0.0, 10000.0, 2, self.setVTop))
        vTopText = QLabel('V<sub>top</sub> (m/s) :')
        self.setVBottom = QLineEdit(
            str(self.dataUI.invData.vBottom), self.groupeOption)
        self.setVBottom.setValidator(
            QtGui.QDoubleValidator(0.0, 10000.0, 2, self.setVBottom))
        vBottomText = QLabel('V<sub>bottom</sub> (m/s) :')
        # self.loadInitModel = QPushButton('Load Initial Model', self.groupeOption)
        minVText = QLabel('Min. velocity (m/s) :')
        self.setVMin = QLineEdit(
            str(self.dataUI.invData.vMin), self.groupeOption)
        self.setVMin.setValidator(
            QtGui.QDoubleValidator(0.0, 10000.0, 2, self.setVMin))
        maxVText = QLabel('Max. velocity (m/s) :')
        self.setVMax = QLineEdit(
            str(self.dataUI.invData.vMax), self.groupeOption)
        self.setVMax.setValidator(
            QtGui.QDoubleValidator(0.0, 10000.0, 2, self.setVMax))
        maxCellText = QLabel('Mesh max. cell size (m²) :')
        self.setMaxCell = QLineEdit(
            str(self.dataUI.invData.meshMaxCellSize), self.groupeOption)
        self.setMaxCell.setValidator(
            QtGui.QDoubleValidator(0.1, 100.0, 2, self.setMaxCell))
        maxDepthText = QLabel('Mesh max. depth (m) :')
        self.setMaxDepth = QLineEdit(
            str(self.dataUI.invData.meshDepthMax), self.groupeOption)
        self.setMaxDepth.setValidator(
            QtGui.QDoubleValidator(5.0, 1000.0, 2, self.setMaxDepth))
        meshQualityText = QLabel('Mesh quality :')
        self.setMeshQuality = QLineEdit(
            str(self.dataUI.invData.meshQuality), self.groupeOption)
        self.setMeshQuality.setValidator(
            QtGui.QDoubleValidator(1.0, 40.0, 2, self.setMeshQuality))
        secNodesText = QLabel('Secondary nodes :')
        self.setSecNodes = QSpinBox(self.groupeOption)
        self.setSecNodes.setRange(0, 12)
        self.setSecNodes.setValue(self.dataUI.invData.secNodes)
        regularizationText = QLabel('Regularization :')
        self.regularizationMode = QComboBox(self.groupeOption)
        self.regularizationMode.addItems(['Smooth (L2)', 'Blocky (L1)'])
        self.refreshMesh = QPushButton('Refresh mesh', self.groupeOption)
        self.refreshMesh.setToolTip(
            'Rebuild and preview the PLC mesh using the current mesh settings.')
        self.runInversion = QPushButton('Run inversion', self.groupeOption)
        groupOptionLayout.addWidget(lambdaText, 0, 0, 1, 1)
        groupOptionLayout.addWidget(self.setLambda, 0, 1, 1, 1)
        groupOptionLayout.addWidget(zWeightText, 0, 2, 1, 1)
        groupOptionLayout.addWidget(self.setZWeight, 0, 3, 1, 1)
        groupOptionLayout.addWidget(vTopText, 1, 0, 1, 1)
        groupOptionLayout.addWidget(self.setVTop, 1, 1, 1, 1)
        groupOptionLayout.addWidget(vBottomText, 1, 2, 1, 1)
        groupOptionLayout.addWidget(self.setVBottom, 1, 3, 1, 1)
        # groupOptionLayout.addWidget(self.loadInitModel, 2, 0, 1, 4)
        # self.loadInitModel.clicked.connect(self._setStartModel)
        groupOptionLayout.addWidget(minVText, 2, 0, 1, 1)
        groupOptionLayout.addWidget(self.setVMin, 2, 1, 1, 1)
        groupOptionLayout.addWidget(maxVText, 2, 2, 1, 1)
        groupOptionLayout.addWidget(self.setVMax, 2, 3, 1, 1)
        groupOptionLayout.addWidget(maxCellText, 3, 0, 1, 1)
        groupOptionLayout.addWidget(self.setMaxCell, 3, 1, 1, 1)
        groupOptionLayout.addWidget(maxDepthText, 3, 2, 1, 1)
        groupOptionLayout.addWidget(self.setMaxDepth, 3, 3, 1, 1)
        groupOptionLayout.addWidget(meshQualityText, 4, 0, 1, 1)
        groupOptionLayout.addWidget(self.setMeshQuality, 4, 1, 1, 1)
        groupOptionLayout.addWidget(secNodesText, 4, 2, 1, 1)
        groupOptionLayout.addWidget(self.setSecNodes, 4, 3, 1, 1)
        groupOptionLayout.addWidget(regularizationText, 5, 0, 1, 1)
        groupOptionLayout.addWidget(self.regularizationMode, 5, 1, 1, 3)
        groupOptionLayout.addWidget(self.refreshMesh, 6, 0, 1, 4)
        groupOptionLayout.addWidget(self.runInversion, 7, 0, 1, 4)
        self.refreshMesh.clicked.connect(self.refreshInversionMesh)
        self.runInversion.clicked.connect(self._runInversion)
        self.groupeOption.setLayout(groupOptionLayout)
        # - Lambda ('lam')
        # - InitialModel
        # - ErrorModel
        # Button for running the inversion
        layout.addWidget(self.groupeOption, 6, 0, 7, 5)
        # self.inversionText = QTextEdit(inversionTab)
        # layout.addWidget(self.inversionText, 10, 0, 1, 5)
        inversionTab.setLayout(layout)
        return inversionTab

    # def _setStartModel(self):
    #     pass

    def _create_inversion_mesh(self):
        """Create the parameter mesh from a PLC using the chosen quality."""
        data = self.dataUI.invData.data
        plc = mt.createParaMeshPLC(
            data, paraDepth=self.dataUI.invData.meshDepthMax,
            paraMaxCellSize=self.dataUI.invData.meshMaxCellSize,
            # The travel-time forward operator does not require an external
            # boundary region. Keeping it creates marker -1 cells after the
            # forward-mesh refinement in recent pyGIMLi releases.
            boundary=0)
        return mt.createMesh(plc, quality=self.dataUI.invData.meshQuality)

    def _read_inversion_mesh_settings(self):
        """Store the current mesh controls before previewing or inverting."""
        self.dataUI.invData.meshMaxCellSize = float(self.setMaxCell.text())
        self.dataUI.invData.meshDepthMax = float(self.setMaxDepth.text())
        self.dataUI.invData.meshQuality = float(self.setMeshQuality.text())
        self.dataUI.invData.secNodes = int(self.setSecNodes.value())
        self.dataUI.invData.blockyModel = (
            self.regularizationMode.currentText() == 'Blocky (L1)')

    def refreshInversionMesh(self):
        """Rebuild and display the mesh without starting an inversion."""
        if self.dataUI.invData.data is None:
            QMessageBox.information(
                self, 'Refresh mesh', 'Load an .sgt file before creating a mesh.')
            return
        try:
            self._read_inversion_mesh_settings()
            self.dataUI.invData.mesh = self._create_inversion_mesh()
        except Exception as exc:
            QMessageBox.warning(self, 'Refresh mesh',
                                f'Unable to create mesh:\n{exc}')
            return
        self.dataUI.meshLoaded = True
        self.dataUI.invData.startModel = None
        self.invModelGraph.axes.clear()
        pgshow(self.dataUI.invData.mesh, ax=self.invModelGraph.axes)
        self.invModelGraph.fig.tight_layout()
        self.invModelGraph.draw()
        self.statusBar.showMessage(
            f'Mesh refreshed: {self.dataUI.invData.mesh.cellCount()} cells, '
            f'quality {self.dataUI.invData.meshQuality:.2f}.', 6000)

    def _runInversion(self):
        """Run a pyGIMLi travel-time inversion from the current UI settings.

        Returns
        -------
        None

        Notes
        -----
        Regularization, velocity limits, mesh settings, and the gradient start
        model are read from the inversion tab. The method rebuilds the mesh,
        runs ``TravelTimeManager.invert``, stores the resulting manager and
        state in ``self.dataUI``, and refreshes both the misfit and model plots.
        """
        # Parameters for inversion:
        self.dataUI.invData.lam = float(self.setLambda.text())
        self.dataUI.invData.zWeight = float(self.setZWeight.text())
        self.dataUI.invData.vTop = float(self.setVTop.text())
        self.dataUI.invData.vBottom = float(self.setVBottom.text())
        self.dataUI.invData.vMin = float(self.setVMin.text())
        self.dataUI.invData.vMax = float(self.setVMax.text())
        # Mesh model and start model:
        self._read_inversion_mesh_settings()
        if self.dataUI.invData.data is not None:
            if self.dataUI.inversionDone:
                self.dataUI.invData.manager = TTMgr(
                    data=self.dataUI.invData.data)
            # Creating mesh with mesh parameters:
            # self.dataUI.invData.data = pg.DataContainer(self.filePicksPath.text())
            # self.dataUI.invData.manager = TTMgr(data=self.dataUI.invData.data)
            # Rebuild even after a preview: this guarantees that the mesh used
            # by the inversion exactly matches the controls at Run time.
            self.dataUI.invData.mesh = self._create_inversion_mesh()
            self.dataUI.invData.startModel = None
            pgshow(self.dataUI.invData.mesh, ax=self.invModelGraph.axes)
            self.invModelGraph.draw()
            # and (self.dataUI.invData.startModel.shape()[0]):
            if (self.dataUI.invData.startModel is None):
                self.dataUI.invData.setStartModelGradient(
                    data=self.dataUI.invData.data, mesh=self.dataUI.invData.mesh)
            # Running the inversion
            self.dataUI.invData.manager.invert(data=self.dataUI.invData.data,
                                               mesh=self.dataUI.invData.mesh,
                                               zWeight=self.dataUI.invData.zWeight,
                                               lam=self.dataUI.invData.lam,
                                               startModel=self.dataUI.invData.startModel,
                                                limits=[
                                                     self.dataUI.invData.vMin, self.dataUI.invData.vMax],
                                                secNodes=self.dataUI.invData.secNodes,
                                                blockyModel=self.dataUI.invData.blockyModel,
                                                verbose=True)
            self.invModelGraph.fig.clear()
            self.invModelGraph.axes = self.invModelGraph.fig.add_subplot(111)
            self.fitGraph.axes.cla()
            drawFirstPicks(ax=self.fitGraph.axes, data=self.dataUI.invData.data, tt=(np.abs(np.asarray(self.dataUI.invData.data(
                't')-np.asarray(self.dataUI.invData.manager.inv.response)))/np.asarray(self.dataUI.invData.data('t')))*100)
            self.fitGraph.axes.set_xlabel('X (m)')
            self.fitGraph.axes.set_ylabel('Data misfit (%)')
            self.fitGraph.fig.tight_layout()
            self.fitGraph.draw()
            _, cBar = self.dataUI.invData.manager.showResult(
                ax=self.invModelGraph.axes, cMap='cividis')
            self.invModelGraph.axes.set_title(
                'Inversion result (chi² = {:.2f}, RMS = {:.6f} s, '
                'RRMS = {:.2f} %)'.format(
                    self.dataUI.invData.manager.inv.chi2(),
                    self.dataUI.invData.manager.inv.absrms(),
                    self.dataUI.invData.manager.inv.relrms()))
            self.invModelGraph.cBar = cBar
            self.dataUI.invData.manager.drawRayPaths(
                ax=self.invModelGraph.axes, color='w', lw=0.3, alpha=0.5)
            self.invModelGraph.fig.tight_layout()
            self.invModelGraph.draw()
            self.dataUI.invData.startModel = None
            self.dataUI.inversionDone = True
            self.dataUI.meshLoaded = False

    def _optional_display_float(self, field):
        """Read an optional numeric display setting, keeping blank as auto."""
        value = field.text().strip()
        if not value:
            return None
        return float(value)

    def _format_display_axes(self, ax, colorbar):
        """Apply user text and font settings after a pyGIMLi plot call."""
        title_size = self.resultTitleFontSize.value()
        label_size = self.resultAxisLabelFontSize.value()
        tick_size = self.resultAxisTickFontSize.value()
        font_family = self.resultFontFamily.currentFont().family()
        # Qt can expose its generic dialog font (e.g. "MS Shell Dlg 2"),
        # which is not an actual Matplotlib font and triggers a warning for
        # every redraw. Fall back to Matplotlib's bundled sans-serif font.
        available_fonts = {font.name for font in font_manager.fontManager.ttflist}
        text_to_render = ''.join(
            [ax.title.get_text(), ax.xaxis.label.get_text(),
             ax.yaxis.label.get_text(), self.resultColorbarTitle.text()])
        if colorbar is not None:
            text_to_render += ''.join(
                label.get_text() for label in
                (colorbar.ax.get_xticklabels() + colorbar.ax.get_yticklabels()))
        try:
            font_path = font_manager.findfont(
                font_manager.FontProperties(family=font_family),
                fallback_to_default=False)
            charmap = font_manager.get_font(font_path).get_charmap()
            has_all_glyphs = all(
                ord(character) in charmap for character in text_to_render
                if character.isprintable())
        except (ValueError, OSError):
            has_all_glyphs = False
        if font_family not in available_fonts or not has_all_glyphs:
            font_family = 'DejaVu Sans'
        title = self.resultTitle.text().strip()
        xlabel = self.resultXAxisTitle.text().strip()
        ylabel = self.resultYAxisTitle.text().strip()
        if title:
            ax.set_title(title)
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)
        ax.title.set_fontsize(title_size)
        ax.xaxis.label.set_fontsize(label_size)
        ax.yaxis.label.set_fontsize(label_size)
        ax.title.set_fontfamily(font_family)
        ax.xaxis.label.set_fontfamily(font_family)
        ax.yaxis.label.set_fontfamily(font_family)
        ax.title.set_fontweight('bold' if self.resultTitleBold.isChecked()
                                else 'normal')
        label_weight = ('bold' if self.resultAxisLabelsBold.isChecked()
                        else 'normal')
        ax.xaxis.label.set_fontweight(label_weight)
        ax.yaxis.label.set_fontweight(label_weight)
        ax.tick_params(labelsize=tick_size)
        for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
            tick_label.set_fontfamily(font_family)
        if colorbar is not None:
            colorbar_title = self.resultColorbarTitle.text().strip()
            colorbar_label = (colorbar.ax.yaxis.label
                              if self.resultColorbarOrientation.currentText() == 'Vertical'
                              else colorbar.ax.xaxis.label)
            label_text = colorbar_title or colorbar_label.get_text()
            colorbar.set_label(
                label_text, fontsize=self.resultColorbarTitleFontSize.value(),
                fontweight=label_weight, fontfamily=font_family)
            # Some pyGIMLi colourbars use an axes title rather than an axis
            # label. Keep it consistent with the editable colourbar label.
            colorbar.ax.title.set_fontsize(
                self.resultColorbarTitleFontSize.value())
            colorbar.ax.title.set_fontfamily(font_family)
            colorbar.ax.title.set_fontweight(label_weight)
            colorbar_label.set_fontfamily(font_family)
            colorbar_label.set_fontweight(label_weight)
            colorbar_label.set_fontsize(self.resultColorbarTitleFontSize.value())
            colorbar.ax.tick_params(labelsize=tick_size)
            for tick_label in (colorbar.ax.get_xticklabels() +
                               colorbar.ax.get_yticklabels()):
                tick_label.set_fontfamily(font_family)

    def _underline_display_title(self, ax):
        """Draw an underline below the title (Matplotlib Text has no underline)."""
        if not self.resultTitleUnderline.isChecked() or not ax.title.get_text():
            return
        canvas = ax.figure.canvas
        canvas.draw()
        renderer = canvas.get_renderer()
        bounds = ax.title.get_window_extent(renderer=renderer)
        transform = ax.transAxes.inverted()
        (x_start, y_line), (x_stop, _) = transform.transform(
            [(bounds.x0, bounds.y0 - 2.0), (bounds.x1, bounds.y0 - 2.0)])
        ax.plot([x_start, x_stop], [y_line, y_line], transform=ax.transAxes,
                clip_on=False, color=ax.title.get_color(),
                linewidth=max(0.6, self.resultTitleFontSize.value()/18.0))

    def updateInversionDisplayControls(self):
        """Disable result-style controls that do not apply to a fit plot."""
        is_fit = self.resultViewSelector.currentText() == 'Observed/modelled fit'
        for widget in (self.resultLogScale, self.resultShowRays,
                       self.resultCoverageMask, self.resultColormap,
                       self.resultColorbarOrientation, self.resultColorMin,
                       self.resultColorMax, self.resultColorbarTitle):
            widget.setEnabled(not is_fit)
        # Rays and coverage masking only make sense for the velocity model.
        is_model = self.resultViewSelector.currentText() == 'Velocity model'
        self.resultShowRays.setEnabled(is_model)
        self.resultCoverageMask.setEnabled(is_model)

    def refreshInversionDisplay(self):
        """Render the selected inversion-result view with the chosen options.

        Returns
        -------
        None

        Notes
        -----
        Depending on the selector, this displays the velocity model,
        standardized coverage, ray coverage, or observed/modelled first picks.
        Colour limits, scale, colormap, colorbar orientation, ray paths, and
        coverage masking are taken from the display controls. Invalid settings
        are reported in the interface without changing the inversion result.
        """
        manager = self.dataUI.invData.manager
        if (not self.dataUI.inversionDone or manager is None or
                not hasattr(manager, 'inv')):
            QMessageBox.information(
                self, 'Display inversion results',
                'Run an inversion before displaying its results.')
            return
        try:
            c_min = self._optional_display_float(self.resultColorMin)
            c_max = self._optional_display_float(self.resultColorMax)
        except ValueError:
            QMessageBox.warning(
                self, 'Display inversion results',
                'Colour limits must be valid numbers, or left blank for automatic limits.')
            return
        if c_min is not None and c_max is not None and c_min >= c_max:
            QMessageBox.warning(
                self, 'Display inversion results',
                'The maximum colour limit must be greater than the minimum.')
            return

        canvas = self.resultDisplayGraph
        canvas.fig.clear()
        ax = canvas.fig.add_subplot(111)
        canvas.axes = ax
        canvas.cBar = None
        view = self.resultViewSelector.currentText()
        colorbar = None
        cmap = self.resultColormap.currentText()
        orientation = self.resultColorbarOrientation.currentText().lower()
        common = {'ax': ax, 'cMap': cmap, 'logScale': self.resultLogScale.isChecked(),
                  'orientation': orientation}
        if c_min is not None:
            common['cMin'] = c_min
        if c_max is not None:
            common['cMax'] = c_max
        try:
            if view == 'Velocity model':
                if self.resultCoverageMask.isChecked():
                    common['coverage'] = manager.standardizedCoverage()
                _, colorbar = manager.showResult(**common)
                if self.resultShowRays.isChecked():
                    manager.drawRayPaths(
                        ax=ax, color='white', lw=0.35, alpha=0.65)
                ax.set_title('Inversion velocity model')
            elif view == 'Standardized coverage':
                coverage = manager.standardizedCoverage()
                common.pop('coverage', None)
                _, colorbar = pgshow(
                    self.dataUI.invData.mesh, coverage, **common)
                ax.set_title('Standardized coverage')
            elif view == 'Ray coverage':
                # Unlike standardized coverage (0/1), ray coverage contains
                # the continuous sensitivity accumulated by the ray paths.
                coverage = manager.rayCoverage()
                common.pop('coverage', None)
                _, colorbar = pgshow(
                    self.dataUI.invData.mesh, coverage, **common)
                ax.set_title('Ray coverage')
            else:
                # pyGIMLi overlays calculated first arrivals on the observed
                # ones when firstPicks=True.
                manager.showFit(firstPicks=True, ax=ax)
                ax.set_title('Observed and modelled first arrivals')
                ax.set_xlabel('X (m)')
                ax.set_ylabel('Traveltime (s)')
        except Exception as exc:
            QMessageBox.warning(
                self, 'Display inversion results',
                f'Unable to render this view:\n{exc}')
            return

        canvas.cBar = colorbar
        self._format_display_axes(ax, colorbar)
        canvas.fig.tight_layout()
        self._underline_display_title(ax)
        canvas.draw()
        self.resultDisplayReady = True

    def saveInversionDisplay(self):
        """Save the currently rendered result figure at publication resolution."""
        if not getattr(self, 'resultDisplayReady', False):
            QMessageBox.information(
                self, 'Save inversion figure',
                'Render a result view before saving it.')
            return
        file_name, _ = QFileDialog.getSaveFileName(
            self, 'Save inversion figure',
            filter=('PNG image (*.png);;JPEG image (*.jpg *.jpeg);;'
                    'TIFF image (*.tif *.tiff);;PDF document (*.pdf);;'
                    'SVG image (*.svg)'))
        if not file_name:
            return
        try:
            self.resultDisplayGraph.fig.savefig(
                file_name, dpi=300, bbox_inches='tight')
            self.statusBar.showMessage(
                f'Inversion figure saved at 300 dpi: {file_name}', 8000)
        except Exception as exc:
            QMessageBox.warning(
                self, 'Save inversion figure', f'Unable to save figure:\n{exc}')

    def _displayInversionResultsTabUI(self):
        """Build the configurable post-inversion result display tab."""
        display_tab = QWidget(self.tabs)
        layout = QGridLayout(display_tab)
        self.resultDisplayGraph = MplCanvas(display_tab, width=9, height=7, dpi=100)
        self.resultDisplayReady = False
        graph_layout = QVBoxLayout()
        graph_layout.addWidget(NavigationToolbar2QT(self.resultDisplayGraph, display_tab))
        graph_layout.addWidget(self.resultDisplayGraph)
        layout.addLayout(graph_layout, 0, 0, 22, 8)

        controls = QGroupBox('Display options', display_tab)
        control_layout = QGridLayout(controls)
        self.resultViewSelector = QComboBox(controls)
        self.resultViewSelector.addItems([
            'Velocity model', 'Standardized coverage', 'Ray coverage',
            'Observed/modelled fit'])
        self.resultColormap = QComboBox(controls)
        self.resultColormap.addItems([
            'cividis', 'viridis', 'plasma', 'plasma_r', 'magma',
            'turbo', 'jet'])
        self.resultColorbarOrientation = QComboBox(controls)
        self.resultColorbarOrientation.addItems(['Vertical', 'Horizontal'])
        self.resultLogScale = QCheckBox('Logarithmic colour scale', controls)
        self.resultShowRays = QCheckBox('Show ray paths', controls)
        self.resultShowRays.setChecked(True)
        self.resultCoverageMask = QCheckBox(
            'Mask model using standardized coverage', controls)
        self.resultColorMin = QLineEdit(controls)
        self.resultColorMax = QLineEdit(controls)
        self.resultColorMin.setPlaceholderText('Auto')
        self.resultColorMax.setPlaceholderText('Auto')
        self.resultTitle = QLineEdit(controls)
        self.resultXAxisTitle = QLineEdit(controls)
        self.resultYAxisTitle = QLineEdit(controls)
        self.resultTitle.setPlaceholderText('Automatic title')
        self.resultXAxisTitle.setPlaceholderText('Automatic x-axis title')
        self.resultYAxisTitle.setPlaceholderText('Automatic y-axis title')
        self.resultColorbarTitle = QLineEdit(controls)
        self.resultColorbarTitle.setPlaceholderText('Automatic colourbar title')
        self.resultTitleFontSize = QSpinBox(controls)
        self.resultAxisLabelFontSize = QSpinBox(controls)
        self.resultAxisTickFontSize = QSpinBox(controls)
        self.resultColorbarTitleFontSize = QSpinBox(controls)
        for selector, value in ((self.resultTitleFontSize, 14),
                                (self.resultAxisLabelFontSize, 11),
                                (self.resultAxisTickFontSize, 9),
                                (self.resultColorbarTitleFontSize, 11)):
            selector.setRange(6, 32)
            selector.setValue(value)
        self.resultFontFamily = QFontComboBox(controls)
        self.resultTitleBold = QCheckBox('Bold title', controls)
        self.resultTitleUnderline = QCheckBox('Underline title', controls)
        self.resultAxisLabelsBold = QCheckBox('Bold axis/colourbar labels', controls)
        self.buttonRefreshInversionDisplay = QPushButton('Apply display options', controls)
        self.buttonSaveInversionDisplay = QPushButton('Save figure (300 dpi)', controls)
        self.buttonRefreshInversionDisplay.clicked.connect(self.refreshInversionDisplay)
        self.buttonSaveInversionDisplay.clicked.connect(self.saveInversionDisplay)
        self.resultViewSelector.currentTextChanged.connect(
            self.updateInversionDisplayControls)

        control_layout.addWidget(QLabel('View:'), 0, 0)
        control_layout.addWidget(self.resultViewSelector, 0, 1)
        control_layout.addWidget(QLabel('Colormap:'), 1, 0)
        control_layout.addWidget(self.resultColormap, 1, 1)
        control_layout.addWidget(QLabel('Colourbar:'), 2, 0)
        control_layout.addWidget(self.resultColorbarOrientation, 2, 1)
        control_layout.addWidget(QLabel('Colour min:'), 3, 0)
        control_layout.addWidget(self.resultColorMin, 3, 1)
        control_layout.addWidget(QLabel('Colour max:'), 4, 0)
        control_layout.addWidget(self.resultColorMax, 4, 1)
        control_layout.addWidget(self.resultLogScale, 5, 0, 1, 2)
        control_layout.addWidget(self.resultShowRays, 6, 0, 1, 2)
        control_layout.addWidget(self.resultCoverageMask, 7, 0, 1, 2)
        control_layout.addWidget(QLabel('Plot title:'), 8, 0)
        control_layout.addWidget(self.resultTitle, 8, 1)
        control_layout.addWidget(QLabel('X-axis title:'), 9, 0)
        control_layout.addWidget(self.resultXAxisTitle, 9, 1)
        control_layout.addWidget(QLabel('Y-axis title:'), 10, 0)
        control_layout.addWidget(self.resultYAxisTitle, 10, 1)
        control_layout.addWidget(QLabel('Colourbar title:'), 11, 0)
        control_layout.addWidget(self.resultColorbarTitle, 11, 1)
        control_layout.addWidget(QLabel('Colourbar-title font size:'), 12, 0)
        control_layout.addWidget(self.resultColorbarTitleFontSize, 12, 1)
        control_layout.addWidget(QLabel('Title font size:'), 13, 0)
        control_layout.addWidget(self.resultTitleFontSize, 13, 1)
        control_layout.addWidget(QLabel('Axis-label font size:'), 14, 0)
        control_layout.addWidget(self.resultAxisLabelFontSize, 14, 1)
        control_layout.addWidget(QLabel('Tick font size:'), 15, 0)
        control_layout.addWidget(self.resultAxisTickFontSize, 15, 1)
        control_layout.addWidget(QLabel('Font family:'), 16, 0)
        control_layout.addWidget(self.resultFontFamily, 16, 1)
        control_layout.addWidget(self.resultTitleBold, 17, 0, 1, 2)
        control_layout.addWidget(self.resultTitleUnderline, 18, 0, 1, 2)
        control_layout.addWidget(self.resultAxisLabelsBold, 19, 0, 1, 2)
        control_layout.addWidget(self.buttonRefreshInversionDisplay, 20, 0, 1, 2)
        control_layout.addWidget(self.buttonSaveInversionDisplay, 21, 0, 1, 2)
        self.updateInversionDisplayControls()
        layout.addWidget(controls, 0, 8, 22, 4)
        return display_tab

    def _modelTabUI(self):
        '''In this tab , we will propose to draw the hodochrones on top
         of the picking and build the corresponding layered model.
        Those models are build using the intercept time-method.
        '''
        modellingTab = QWidget()
        # Pick loading
        self.filePicksPath = QLabel('File path', modellingTab)
        self.filePicksPath.setAlignment(QtCore.Qt.AlignCenter)
        self.buttonGetSgtFile = QPushButton('...', modellingTab)
        self.buttonGetSgtFile.clicked.connect(self._loadPicking)
        # Graph with the hodochrones:
        self.hodochronesGraph = MplCanvas(modellingTab)
        hodochronesToolbar = NavigationToolbar2QT(
            self.hodochronesGraph, modellingTab)
        hodochronesWidget = QVBoxLayout()
        hodochronesWidget.addWidget(
            hodochronesToolbar, alignment=QtCore.Qt.AlignLeft)
        hodochronesWidget.addWidget(self.hodochronesGraph)
        # Graph with the model
        self.modelGraph = MplCanvas(modellingTab)
        modelToolbar = NavigationToolbar2QT(self.modelGraph, modellingTab)
        modelWidget = QVBoxLayout()
        modelWidget.addWidget(modelToolbar, alignment=QtCore.Qt.AlignRight)
        modelWidget.addWidget(self.modelGraph)
        # Options for the tab:
        self.groupOptionModelling = QGroupBox()
        self.groupOptionModelling.setTitle('Options')
        self.sourceSelector = QComboBox(self.groupOptionModelling)
        self.receiversSelector = QComboBox(self.groupOptionModelling)
        labelSourceReceiver = QLabel(
            'Select source and orientation :', self.groupOptionModelling)
        layoutOptions = QGridLayout()
        layoutOptions.addWidget(labelSourceReceiver, 1,
                                0, 1, 4, alignment=QtCore.Qt.AlignRight)
        layoutOptions.addWidget(self.sourceSelector, 1, 4, 1, 4)
        layoutOptions.addWidget(self.receiversSelector, 1, 8, 1, 2)
        self.movePoints = QPushButton('Move Points', self.groupOptionModelling)
        self.movePoints.setCheckable(True)
        self.movePoints.clicked.connect(self.setAnimationModelling)
        self.nbLayersSelector = QSpinBox(self.groupOptionModelling)
        self.nbLayersSelector.setMaximum(3)
        self.nbLayersSelector.setMinimum(2)
        self.nbLayersSelector.setValue(self.dataUI.modellingData.nbLayers)
        self.nbLayersSelector.valueChanged.connect(self._initModelling)
        nbLayersLabel = QLabel('Select the number of layers :')
        layoutOptions.addWidget(nbLayersLabel, 0, 0, 1,
                                5, alignment=QtCore.Qt.AlignRight)
        layoutOptions.addWidget(self.nbLayersSelector, 0, 5, 1, 5)
        layoutOptions.addWidget(self.movePoints, 2, 0, 1, 10)
        self.groupOptionModelling.setLayout(layoutOptions)
        # Setup of the layout:
        layout = QGridLayout(modellingTab)
        layout.addWidget(self.filePicksPath, 0, 0, 1, 9)
        layout.addWidget(self.buttonGetSgtFile, 0, 9, 1, 1)
        layout.addLayout(hodochronesWidget, 1, 0, 7, 5)
        layout.addLayout(modelWidget, 1, 5, 7, 5)
        layout.addWidget(self.groupOptionModelling, 9, 0, 2, 5)
        modellingTab.setLayout(layout)
        self.aniModelling = None
        return modellingTab

    def saveStatePicking(self):
        """Serialize the current data and picking session to a pickle file."""
        fName, _ = QFileDialog.getSaveFileName(
            self, 'Select file to save', filter='Pickled structure (*.pkl)')
        if fName != "":
            dataLoaded = picklingStatus()
            dataLoaded.paths = self.dataUI.paths
            dataLoaded.geometry = self.dataUI.geometry
            dataLoaded.sisData = self.dataUI.sisData
            dataLoaded.beginTime = self.dataUI.beginTime
            dataLoaded.picking = self.dataUI.picking
            dataLoaded.pickingError = self.dataUI.pickingError
            dataLoaded.manualPickingMask = self.dataUI.manualPickingMask
            dataLoaded.manualDefectiveOverrides = (
                self.dataUI.manualDefectiveOverrides)
            dataLoaded.sisFileId = self.dataUI.sisFileId
            file = open(fName, 'wb')
            pickle.dump(dataLoaded, file)
            file.close()
            self.statusBar.showMessage(f'Result saved to {fName}', 10000)
        else:
            self.statusBar.showMessage('Result was NOT saved!', 2000)

    def loadStatePicking(self):
        """Restore a previously serialized picking session into the interface."""
        fName, _ = QFileDialog.getOpenFileName(
            self, 'Select file to load', filter='Pickled structure (*.pkl)')
        if fName != "":
            file = open(fName, 'rb')
            dataLoaded = pickle.load(file)
            file.close()
            # Assign the different elements :
            if len(dataLoaded.paths.directory) != 0:
                if self.dataUI.dataLoaded:
                    self.aniZoom.event_source.stop()
                    self.aniMain.event_source.stop()
                    self.mainGraph.mpl_disconnect(self.connectMouse)
                    self.mainGraph.mpl_disconnect(self.connectPress)
                    self.mainGraph.mpl_disconnect(self.connectRelease)
                self.dataUI.paths = dataLoaded.paths
                self.dataUI.geometry = dataLoaded.geometry
                self.dataUI.sisData = dataLoaded.sisData
                self.dataUI.beginTime = dataLoaded.beginTime
                self.dataUI.picking = dataLoaded.picking
                self.dataUI.pickingError = dataLoaded.pickingError
                self.dataUI.manualPickingMask = getattr(
                    dataLoaded, 'manualPickingMask',
                    np.isfinite(dataLoaded.picking).astype(bool))
                self.dataUI.manualDefectiveOverrides = getattr(
                    dataLoaded, 'manualDefectiveOverrides',
                    np.zeros_like(dataLoaded.picking, dtype=np.int8))
                self.dataUI.sisFileId = dataLoaded.sisFileId
                self.dataUI.dataLoaded = True
                # Change the interface:
                self.spinBoxCurrSelect.setMinimum(0)
                self.spinBoxCurrSelect.setMaximum(
                    len(self.dataUI.sisData[0])-1)
                self.spinBoxCurrSelect.setPrefix('Trace number ')
                self.spinBoxCurrSelect.valueChanged.connect(
                    self.traceNumberChanged)
                self.spinBoxCurrSelect.setEnabled(True)
                self.buttonManageDefective.setEnabled(True)
                self.dataUI.animationPicking.changedSelect = True
                self.updateTab0()
            else:
                self.statusBar.showMessage('Empty state loaded', 2000)
        else:
            self.statusBar.showMessage('No status loaded!', 2000)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    window = Window()
    window.show()
    sys.exit(app.exec_())
