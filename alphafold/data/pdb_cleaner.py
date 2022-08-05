from absl import flags
from absl import app
import numpy as np
import string

flags.DEFINE_string('pdb_path_input', None, 'Path to the pdb file to clean.')
flags.DEFINE_string('pdb_path_output', None, 'Path to the pdb file to the cleaned.')
FLAGS = flags.FLAGS

def main(argv):
    with open(FLAGS.pdb_path_input, "r") as f:
        handle = f.readlines()

    unique_identifier = []
    nb_TER = 0
    for line in handle:
        if "ATOM" in line:
            if line[20] != ' ':
                line_list = list(line)
                unique_identifier.append(line_list[20])
            else:
                line_list = list(line)
                unique_identifier.append(line_list[21])

        if line.strip("\n") == "TER":
            nb_TER += 1

    chain_number = 0
    all_lines = []
    for line in handle:
        if "ATOM" in line:
            line_list = list(line)
            line_list[21] = string.ascii_uppercase[chain_number]
            line_list[20] = ' '
            line = ''.join(line_list)

        elif line.strip("\n") == "TER":
            chain_number += 1


        all_lines.append(line)

    with open(FLAGS.pdb_path_output, "w") as f:
        f.writelines(all_lines)

if __name__ == '__main__':
  flags.mark_flags_as_required([
      'pdb_path_input',
      'pdb_path_output'
  ])

  app.run(main)