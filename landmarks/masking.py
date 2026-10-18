"""Independent image-only Otsu baseline. No coordinates are changed."""
import numpy as np
from scipy.ndimage import label
from skimage.filters import threshold_otsu


def otsu_mask(image, valid, tissue_is_bright, max_speck_pixels=25):
    image=np.asarray(image)
    valid=np.asarray(valid)
    if image.ndim!=2 or valid.shape!=image.shape or valid.dtype!=bool:
        raise ValueError('Expected 2D image and matching Boolean support mask')
    values=image[valid]
    if not len(values) or not np.isfinite(values).all() or np.ptp(values)==0:
        raise ValueError('Valid support must contain finite, nonconstant intensities')
    if max_speck_pixels<0: raise ValueError('Speck limit must be nonnegative')
    threshold=float(threshold_otsu(values,nbins=256))
    def classify(t):
        return valid & ((image>t) if tissue_is_bright else (image<t))
    raw=classify(threshold)
    # Eight-connected components; no hole filling, closing, dilation or erosion.
    components,count=label(raw,structure=np.ones((3,3),int))
    counts=np.bincount(components.ravel())
    remove=(counts<=max_speck_pixels); remove[0]=False
    cleaned=raw & ~remove[components]
    delta=.02*float(np.ptp(values))
    low,high=classify(threshold-delta),classify(threshold+delta)
    stats=dict(method='Otsu, 256 histogram bins',threshold=threshold,
        tissue_rule='intensity > threshold' if tissue_is_bright else 'intensity < threshold',
        valid_pixels=int(valid.sum()),raw_tissue_pixels=int(raw.sum()),
        cleaned_tissue_pixels=int(cleaned.sum()),removed_pixels=int((raw & ~cleaned).sum()),
        components_before=int(count),max_removed_component_pixels=max_speck_pixels,
        sensitivity_delta=.02,sensitivity_delta_definition='2% of valid intensity range, not anatomical uncertainty',
        sensitivity_delta_intensity=delta,
        pixels_changed_at_lower_threshold=int(np.count_nonzero(low!=raw)),
        pixels_changed_at_upper_threshold=int(np.count_nonzero(high!=raw)))
    return dict(raw=raw,cleaned=cleaned,lower_threshold=low,upper_threshold=high),stats


def histology_hysteresis(image, valid, weak_fraction=.005, strong_fraction=.05):
    """Grow dark tissue seeds into pale connected tissue; never fill holes.

    Background = median of brightest valid decile. Threshold margins are
    relative to valid intensity range, with a robust background-noise floor.
    Fixed starting parameters, not fitted to reference-mask overlap.
    """
    from scipy.ndimage import binary_propagation
    image=np.asarray(image); valid=np.asarray(valid)
    if image.ndim!=2 or valid.shape!=image.shape or valid.dtype!=bool:
        raise ValueError('Expected 2D image and matching Boolean support mask')
    values=image[valid]
    if not len(values) or not np.isfinite(values).all() or np.ptp(values)==0:
        raise ValueError('Need finite, nonconstant valid intensities')
    if not 0<weak_fraction<strong_fraction<1:
        raise ValueError('Require 0 < weak fraction < strong fraction < 1')
    span=float(np.ptp(values))
    brightest=values[values>=np.percentile(values,90)]
    background=float(np.median(brightest))
    noise=float(1.4826*np.median(np.abs(brightest-background)))
    weak_margin=max(weak_fraction*span,3*noise)
    strong_margin=max(strong_fraction*span,2*weak_margin)
    darkness=background-image
    possible=valid & (darkness>weak_margin)
    strong=valid & (darkness>strong_margin)
    if not strong.any(): raise ValueError('No strong tissue seeds; inspect background estimate')
    result=binary_propagation(strong,mask=possible,structure=np.ones((3,3),bool))
    stats=dict(method='Background-relative hysteresis, 8-connected propagation',
        background_estimator='Median of brightest 10% of valid pixels',
        estimated_background=background,bright_tail_robust_sigma=noise,
        intensity_span=span,weak_fraction=weak_fraction,strong_fraction=strong_fraction,
        weak_darkness_margin=weak_margin,strong_darkness_margin=strong_margin,
        possible_tissue_intensity_below=background-weak_margin,
        strong_tissue_intensity_below=background-strong_margin,
        strong_seed_pixels=int(strong.sum()),possible_tissue_pixels=int(possible.sum()),
        hysteresis_tissue_pixels=int(result.sum()),pale_connected_pixels_added=int((result&~strong).sum()),
        disconnected_possible_pixels_rejected=int((possible&~result).sum()),
        assumptions='Near-white background occupies the bright tail; uniform background in this processed BigBrain slice',
        limitations='Fixed initial margins require validation; pale disconnected tissue may be missed; noise floor estimated from truncated bright tail is not a calibrated uncertainty')
    return dict(strong=strong,possible=possible,mask=result),stats
