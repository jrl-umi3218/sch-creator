import bpy
import sys

argv = sys.argv
argv = argv[argv.index("--") + 1:] # get all args after "--"

fbx_in = argv[0]
fbx_out = argv[1]

bpy.ops.import_scene.obj(filepath=fbx_in, axis_forward='Y', axis_up='Z')
bpy.ops.export_scene.qc(filepath=fbx_out)

