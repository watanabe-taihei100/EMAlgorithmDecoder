import sys
sys.path.append("../main")
import EMdecoderPolarization

num_patterns = 1
attenuation = 0.0

decode_runner = EMdecoderPolarization.DecodeRunnerWithPolarization(num_patterns)

pattern_id = 0
pattern_file_name = "pattern_0.txt"
decode_runner.read_pattern_file(pattern_file_name, pattern_id)
decode_runner.read_encoded_image_txt("./encoded_image_Htype.txt", pattern_id, 0)
decode_runner.read_encoded_image_txt("./encoded_image_Vtype.txt", pattern_id, 90)

start_index = 0

outdirname = "test2"
decode_runner.construct_matrix(attenuation)
decode_runner.run(outdirname, 10, start_index=start_index)