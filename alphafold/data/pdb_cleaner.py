from absl import flags
from absl import app

flags.DEFINE_string('pdb_path_input', None, 'Path to the pdb file to clean.')
flags.DEFINE_string('pdb_path_output', None, 'Path to the pdb file to clean.')
FLAGS = flags.FLAGS

def main(argv):
    with open(FLAGS.pdb_path_input, "r") as f:
        handle = f.readlines()

    all_lines = []
    for line in handle:
        if "ATOM" in line:
            if line[20] != ' ':
                line_list = list(line)
                line_list[21] = line_list[20]
                line_list[20] = ' '
                line = ''.join(line_list)

        all_lines.append(line)

    with open(FLAGS.pdb_path_output, "w") as f:
        f.writelines(all_lines)

if __name__ == '__main__':
  flags.mark_flags_as_required([
      'pdb_path_input',
      'pdb_path_output'
  ])

  app.run(main)