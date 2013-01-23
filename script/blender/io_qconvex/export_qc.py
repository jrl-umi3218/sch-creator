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

__author__ = ["Joris Vaillant"]
__version__ = '0.1'
__bpydoc__ = """\
This script exports a Mesh to a qc qconvex cloud file format.

"""

import bpy

def save(filepath):
  cloud = []
  for obj in bpy.context.scene.objects:
    if obj.type == 'MESH':
      matrix = obj.matrix_world
      cloud.extend([vert.co * matrix for vert in obj.data.vertices])

  # write the faces to a file
  with open(filepath, 'w') as file:
    file.write('3\n')
    file.write('%s\n' % len(cloud))
    data = ['%s %s %s\n' % (c.x, c.y, c.z) for c in cloud]
    file.write(''.join(data))

  return {'FINISHED'}

