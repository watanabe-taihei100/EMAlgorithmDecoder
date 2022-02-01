from ..main import EMdecoderPolarization
import json

num_patterns = 8
attenuation = 0.02
dirname = "scan_16keV_blend_33"
seed = 0


def read_corr_factor():
    filename = "/Users/watanabe/work/SPring8_analysis/polarization_results/16keV_result.json"
    ret = [0, 0]
    with open(filename, "r") as f:
        obj = json.load(f)
        ret[0] = obj["correction_factor"]["corr_Htype"]
        ret[1] = obj["correction_factor"]["corr_Vtype"]
    return ret

decode_runner = EMdecoderPolarization.DecodeRunnerWithPolarization(num_patterns)

# calculate weight
E_count = {}
# min_count = 10**18
for pattern_id in range(num_patterns):
    for theta_z in ["0deg", "90deg"]:
        event_count = 0
        for eventtype in ["Htype", "Vtype"]:
            filename = "/Users/watanabe/work/SPring8_analysis/products/imaging/16keV/{}/scan_all/encoded_image_{}/pattern_{}.txt".format(theta_z, eventtype, pattern_id)
            with open(filename, "r") as f:
                for row in f:
                    event_count += sum(map(int, row.rstrip().split()))
            E_count["{}_{}".format(theta_z, pattern_id)] = event_count

E_count_all = sum(E_count.values())
for key in E_count.keys():
    E_count[key] /= (E_count_all / len(E_count))

print(E_count)
# print(min_count)
corr_factor = read_corr_factor()
min_count = 0


for pattern_id in range(num_patterns):
    pattern_file_name = "/Users/watanabe/work/Coded_aperture/patterns/pattern_{}_2D.txt".format(pattern_id)
    decode_runner.read_pattern_file(pattern_file_name, pattern_id)

    encoded_image_filename_bases = [
            "/Users/watanabe/work/SPring8_analysis/products/imaging/16keV/0deg/scan_all/encoded_image_xm_yp_frac_5",
            "/Users/watanabe/work/SPring8_analysis/products/imaging/16keV/90deg/scan_all/encoded_image_xm_yp_frac_5",
            "/Users/watanabe/work/SPring8_analysis/products/imaging/16keV/90deg/scan_all/encoded_image_xp_yp",
            "/Users/watanabe/work/SPring8_analysis/products/imaging/16keV/0deg/scan_all/encoded_image_xm_ym",
            "/Users/watanabe/work/SPring8_analysis/products/imaging/16keV/0deg/scan_all/encoded_image_xp_ym_frac_7",
            "/Users/watanabe/work/SPring8_analysis/products/imaging/16keV/90deg/scan_all/encoded_image_xp_ym_frac_3"
    ]

    for filename_base in encoded_image_filename_bases:
        theta_z = "90deg" if "/90deg/" in filename_base else "0deg"
        beam_weight = E_count["{}_{}".format(theta_z, pattern_id)]
        weight_Htype = beam_weight / corr_factor[0]
        weight_Vtype = beam_weight / corr_factor[1]
        print(weight_Htype, weight_Vtype)
        decode_runner.read_encoded_image_txt_double(filename_base, pattern_id, weight_Htype=weight_Htype, weight_Vtype=weight_Vtype, sample=min_count, seed=seed)

start_index = 0

decode_runner.construct_matrix(attenuation)
decode_runner.run(dirname, 1500, start_index=start_index)