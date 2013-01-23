import bpy
import sys

sce = bpy.context.scene

argv = sys.argv
argv = argv[argv.index("--") + 1:] # get all args after "--"

fbx_in = argv[0]
fbx_out = argv[1]

bpy.ops.import_scene.x3d(filepath=fbx_in, axis_forward='Y', axis_up='Z')

for obj in bpy.context.scene.objects:
  if obj.type == 'MESH':
    sce.objects.active = obj
    obj.modifiers.new(name='toto', type='DECIMATE')
    obj.modifiers['toto'].ratio = 0.01
    bpy.ops.object.modifier_apply(modifier='toto')

bpy.ops.export_scene.qc(filepath=fbx_out)

