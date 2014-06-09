import os
import sys

import bpy


# file loader per extension
loader_factory = {
  '.obj': lambda file: bpy.ops.import_scene.obj(filepath=file, axis_forward='Y', axis_up='Z'),
  '.wrl': lambda file: bpy.ops.import_scene.x3d(filepath=file, axis_forward='Y', axis_up='Z')
}


def convert(input_filename, output_filename):
  "convert a file to a qc file."
  ext = os.path.splitext(input_filename)[1]
  try:
    loader = loader_factory[ext]
  except KeyError:
    raise RuntimeError('%s is not supported' % ext)

  loader(input_filename)
  bpy.ops.export_scene.qc(filepath=output_filename)


def main():
  argv = sys.argv
  argv = argv[argv.index("--") + 1:] # get all args after "--"

  fbx_in = argv[0]
  fbx_out = argv[1]
  convert(fbx_in, fbx_out)


if __name__ == '__main__':
  main()
