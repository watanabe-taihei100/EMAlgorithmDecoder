import sys
sys.path.append("../main")
import EMdecoderSimple

num_patterns = 1
attenuation = 0.0

decode_runner = EMdecoderSimple.DecodeRunner(num_patterns)

pattern_id = 0
pattern_file_name = "pattern_0.txt"
encoded_image_filename = "encoded_image.txt"
decode_runner.read_pattern_file(pattern_file_name, pattern_id)
decode_runner.read_encoded_image_txt(encoded_image_filename, pattern_id)

outdirname = "test1"
decode_runner.construct_matrix(attenuation)
beta = 1
decode_runner.run(outdirname, 5, beta, start_index=0)


