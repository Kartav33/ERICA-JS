from scipy.ndimage import grey_opening
from scipy.stats import linregress
import numpy as np
import pandas as pd

#-----------------------------------------------------------------------
# To calculate the primary baseline
#-----------------------------------------------------------------------

def lulu_baseline(intensity, burst_peak, window_hr, smooth_hr, cadence):

    signal = intensity.loc[burst_peak - pd.Timedelta(hours=window_hr):burst_peak + pd.Timedelta(hours=window_hr)].copy()
    n = max(1, round(pd.Timedelta(hours=smooth_hr) / pd.Timedelta(cadence)))

    # Ensure regular cadence so missing timestamps are also detected
    signal = signal.asfreq(cadence)

    # Fill only short internal gaps (up to 3 samples)
    filled = signal.interpolate(method="time", limit=3, limit_area="inside")
    baseline = pd.Series(np.nan, index=filled.index)

    # Identify continuous segments without missing values
    valid = filled.notna()
    groups = valid.ne(valid.shift()).cumsum()

    for _, segment in filled.groupby(groups):

        if segment.isna().any():
            continue

        if len(segment) < n + 1:
            continue

        baseline.loc[segment.index] = grey_opening(segment.to_numpy(), size=n + 1)

    return baseline

#-----------------------------------------------------------------------
# Shifting the primary baseline by an offset rule to find the threshold
#-----------------------------------------------------------------------

def shift_baseline(intensity, baseline_fit, reference_time, background_hr, quantile, direction="before"):

    signal = intensity.reindex(baseline_fit.index)
    window = pd.Timedelta(hours=background_hr)

    if direction == "before":
        background = signal.loc[reference_time - window:reference_time]
    elif direction == "after":
        background = signal.loc[reference_time:reference_time + window]
    else:
        raise ValueError("direction must be 'before' or 'after'")

    difference = (background - baseline_fit.reindex(background.index)).dropna()

    if difference.empty:
        return baseline_fit * np.nan, np.nan

    offset = difference.quantile(quantile)

    return baseline_fit + offset, offset

#-----------------------------------------------------------------------
# Detecting the ending point, end_flags = 0 and 1
#-----------------------------------------------------------------------

def detect(threshold, signal, burst_peak, burst_base, cadence_td, rotation):

    for t in signal.loc[burst_peak:].index:
        if pd.isna(signal.loc[t]) or signal.loc[t] > threshold.loc[t]:
            continue

        duration = t - burst_base
        flag = int(duration < rotation)

        confirmation = rotation * (1 + flag)
        check_end = t + confirmation

        if check_end > signal.index[-1]:
            continue

        window = signal.loc[t:check_end]
        threshold_window = threshold.loc[t:check_end]

        expected = int(confirmation / cadence_td) + 1

        if len(window) < expected or window.isna().any():
            continue

        if (window > threshold_window).any():
            continue

        return t, flag

    return pd.NaT, 0

#-----------------------------------------------------------------------
# Fallback case, in case the code never finds a natural ending. end_flag = 2
#-----------------------------------------------------------------------

def detect_fallback(threshold, signal, burst_peak, cadence_td):

    longest_duration = pd.Timedelta(0)
    active_end = pd.NaT

    for t in signal.loc[burst_peak:].index:
        if pd.isna(signal.loc[t]) or signal.loc[t] > threshold.loc[t]:
            #print(" I am here")
            continue

        # Locating the longest continuous quiet interval
        win = signal.loc[t:]
        quiet = win.notna() & (win <=threshold.reindex(win.index))

        # Stopping at the onset of the identified quite interval
        end_idx =  np.flatnonzero(~quiet.to_numpy())
        n_quiet = end_idx[0] if len(end_idx) else len(quiet)

        duration = n_quiet * cadence_td
        if duration > longest_duration:
            longest_duration = duration
            active_end = t

    return active_end, 2

#-----------------------------------------------------------------------
# Finding the end point, nested function
#-----------------------------------------------------------------------

def find_active_end_linear(intensity, baseline_fit, burst_peak, burst_base, background_hr, quantile, rotation_hr, cadence):

    signal = intensity.reindex(baseline_fit.index)

    rotation = pd.Timedelta(hours=rotation_hr)
    cadence_td = pd.Timedelta(cadence)

    # CASE 0: Pre-burst background
    threshold_used, offset = shift_baseline(intensity, baseline_fit, burst_base, background_hr, quantile, direction="before")
    active_end, flag = detect(threshold_used, signal, burst_peak, burst_base, cadence_td, rotation)
    threshold_flag = 0

    # CASE 1: Post-peak background
    if pd.isna(active_end):

        threshold_used, offset = shift_baseline(intensity, baseline_fit, burst_peak, background_hr, quantile, direction="after")
        active_end, flag = detect(threshold_used, signal, burst_peak, burst_base, cadence_td, rotation)
        threshold_flag = 1

    # CASE 2: Longest quiet interval
    if pd.isna(active_end):

        active_end, flag = detect_fallback(threshold_used, signal, burst_peak, cadence_td)
        threshold_flag = 2

    # Looking for the point closest to the baseline to find the true end of the signal
    window = signal.loc[active_end:active_end +rotation]
    difference = (window - baseline_fit.reindex(window.index)).dropna()
    if not difference.empty:
        active_end = difference.idxmin()

    return active_end, flag, threshold_flag, threshold_used

#-----------------------------------------------------------------------
# Detecting the the active start based on the same threshold used to detect the ending point
#-----------------------------------------------------------------------

def find_active_start_linear(intensity, hourly_intensity, threshold_fit, prominent_base,
                             lookback_hr, rotation_hr, slope, cadence
                             ):

    rotation = pd.Timedelta(hours=rotation_hr)
    cadence_td = pd.Timedelta(cadence)

    # Extend the existing shifted threshold backward
    def extend_threshold(index):
        t = (index - threshold_fit.index[0]).total_seconds() / 3600
        return pd.Series(
            threshold_fit.iloc[0] + slope * t,
            index=index
        )

    # Original backward detector, reusable for cases 1 and 2
    def backward_detect(anchor):

        index = pd.date_range(
            end=anchor,
            periods=int(np.ceil(
                pd.Timedelta(hours=lookback_hr) / cadence_td
            )) + 1,
            freq=cadence
        )

        signal = intensity.reindex(index)
        threshold = extend_threshold(index)
        difference = threshold - signal

        # Half a Jovian rotation of continuous quiet emission
        n = int(np.ceil((rotation / 2) / cadence_td))
        count = 0

        provisional = pd.NaT

        for t, value in difference.iloc[::-1].items():

            if pd.notna(value) and value >= 0:
                count += 1
            else:
                count = 0

            if count >= n:
                active = difference.loc[t:anchor]
                active = active[active < 0]

                if not active.empty:
                    provisional = active.index[0]

                break

        if pd.isna(provisional):
            return pd.NaT, pd.NaT

        # Minimum-residual refinement, one rotation backward
        window = signal.loc[provisional - rotation:provisional]
        residual = (
            window - threshold.reindex(window.index)
        ).dropna()

        refined = residual.idxmin() if not residual.empty else provisional

        return provisional, refined

    # ============================================================
    # CASE 1: Standard backward detection
    # ============================================================

    provisional_start_1, active_start_1 = backward_detect(prominent_base)

    if pd.isna(active_start_1):
        return prominent_base, 0, pd.NaT

    active_start = active_start_1
    active_start_2 = pd.NaT
    flag = 1

    # ============================================================
    # CASE 2: Inspect earlier activity using hourly derivative
    # ============================================================

    # W2: anchored to provisional start, ending at refined start
    w2_start = provisional_start_1 - 2 * rotation
    w2_end = active_start_1

    hourly = hourly_intensity.loc[w2_start:w2_end].dropna()

    if len(hourly) >= 3:

        # Derivative in chronological time
        derivative = hourly.diff()

        # Local maxima: derivative changes from positive to negative
        peaks = hourly.index[
            (derivative.shift(1) > 0) & (derivative < 0)
        ]

        for peak_time in peaks:

            # Peak itself is the sample immediately before the sign change
            peak_pos = hourly.index.get_loc(peak_time)
            if peak_pos == 0:
                continue

            earlier_peak = hourly.index[peak_pos - 1]

            # Require two consecutive hourly samples above threshold
            peak_window = hourly.loc[
                earlier_peak:earlier_peak + pd.Timedelta(hours=1)
            ]

            threshold_window = extend_threshold(peak_window.index)

            if len(peak_window) < 2:
                continue

            if not (peak_window > threshold_window).all():
                continue

            # Restart original backward detector from earlier peak
            provisional_start_2, candidate_start_2 = backward_detect(
                earlier_peak
            )

            if (
                pd.notna(candidate_start_2)
                and candidate_start_2 < active_start_1
            ):
                active_start_2 = candidate_start_2
                #active_start = active_start_2
                flag = 2
                break

    return active_start, flag, active_start_2
