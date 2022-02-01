import sys
sys.path.append("../main")
import EMdecoderSimple

num_patterns = 8
attenuation = 0.02
dirname = "scansim_16keV_0deg_center_double"

decode_runner = EMdecoderSimple.DecodeRunner(num_patterns)

for pattern_id in range(num_patterns):
    encoded_image_filename = "/Users/watanabe/work/SPring8_simulation/scan/products/16keV/0deg/random/pattern{}/encoded_image_center.txt".format(pattern_id)

    pattern_file_name = "/Users/watanabe/work/Coded_aperture/patterns/pattern_{}_2D.txt".format(pattern_id)
    decode_runner.read_pattern_file(pattern_file_name, pattern_id)
    decode_runner.read_encoded_image_txt(encoded_image_filename, pattern_id)

decode_runner.construct_matrix(attenuation)
beta = 1
decode_runner.run(dirname, 1000, beta, start_index=0)


