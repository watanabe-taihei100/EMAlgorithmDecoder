# Decode Runner using EM Algorithm for Coded Aperture Imaging

## Scripts

- main

||description|
|---|---|
|EMdecoderSimple.py|simple imaging|
|EMdecoderPolarization.py|imaging with polarization|

- example

||description|
|---|---|
|run_simple.py|example to run EMdecoderSimple.py|
|run_plarization.py|example to run EMdecoderPolarization.py|

- cpp

used for drawing image with ROOT.


## Main Sequence

1. Set parameters

Use `read_pattern_file()` and `read_encoded_image_txt()`.

2. Construct matrix

Use `construct_matrix`.

3. Run

Use `run()`.
`E_step()` and `M_step()` are main algorithm.


## Parameters and Variables

### in EMdecoderSimple.py

parameters need to set is as follows.

||description|example in SPring8_2021B experiment|
|---|---|---|
|num_patterns|the number of patterns|8|
|Ndx,Ndy|the number of detector pixels per pattern|896|
|Nmx,Nmy|the number of pattern elements per pattern|64|
|Nsx,Nsy|the number of sky elements |51|
|mask_pitch|length between pattern elements (um)|35|
|pixel_pitch|length between detector pixels (um)|2.5|
|sky_pitch|length between sky elements (arcsec)|20|
|detector_to_mask|distance from detector to mask (um)|20|

variables used in EM Algorithm is as follows.

||description|
|---|---|
|Mat|matrix from sky to detector|
|Pattern|mask patterns|
|EncodedImage|events detected by experiment|
|detector_X, detector_Y|detector coordinate|
|sky_X, sky_Y|sky coordinate|
|S_dist|expected sky destribution|
|D_dist|expected detector destribution|


### additional parameters in EMdecoderPolarization.py

||description|example in SPring8_2021B experiment|
|---|---|---|
|Ndp|the number of angle|2|
|angle_pitch|pitch of angle (degree)|90|
|modulation_factor|modulation factor of detector|0.103..|
