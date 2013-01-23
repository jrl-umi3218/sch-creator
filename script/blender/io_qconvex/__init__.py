# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8 compliant>

bl_info = {
  "name": "qconvex cloud format",
  "author": "Joris Vaillant",
  "blender": (2, 5, 7),
  "api": 35622,
  "location": "File > Import-Export",
  "description": "Export to qconvex cloud format",
  "warning": "",
  "support": 'OFFICIAL',
  "category": "Import-Export"}


if "bpy" in locals():
  import imp
  if 'export_qc' in locals():
    imp.reload(export_qc)

import bpy
from bpy.props import StringProperty
from bpy_extras.io_utils import ExportHelper


class ExportQConvex(bpy.types.Operator, ExportHelper):
  '''Save QConvex cloud'''
  bl_idname = "export_scene.qc"
  bl_label = "Export QC"

  filename_ext = '.qc'
  filepath = StringProperty(name="File Path", description="Filepath used for exporting the QC file", maxlen= 1024, default= "", subtype='FILE_PATH')

  def execute(self, context):
    from . import export_qc
    return export_qc.save(self.filepath)

def menu_func_export(self, context):
  self.layout.operator(ExportQConvex.bl_idname, text="QConvex cloud (.qc)")


def register():
  bpy.utils.register_module(__name__)

  bpy.types.INFO_MT_file_export.append(menu_func_export)


def unregister():
  bpy.utils.unregister_module(__name__)

  bpy.types.INFO_MT_file_export.remove(menu_func_export)

# NOTES
# - blender version is hardcoded

if __name__ == "__main__":
  register()

