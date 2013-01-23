#! /usr/bin/python
class Vrml(object):
  def __init__(self):
    self.points = []
    self._isPoint = False

  def load(self, file):
    for l in file.readlines():
      l = l.strip()
      if self._isPoint:
        self._readPoint(l)
      else:
        if 'point' in l:
          self._isPoint = True

  def _readPoint(self, l):
    if ']' in l:
      self._isPoint = False
    else:
      p = l.split(' ')
      p[2] = p[2][0:-1]
      self.points.append(p)
  
class Cloud(object):
  def __init__(self, points):
    self._points = points

  def save(self, file):
    # File format :
    # dimension \n number of points \n points coordinates  
    file.write('3 \n')
    file.write('%s\n' % len(self._points))
    data = ['%s %s %s\n' % (p[0], p[1], p[2]) for p in self._points]
    file.write(''.join(data))

if __name__ == '__main__':
  import sys

  if len(sys.argv) != 3:
    print 'Usage : %s input.wrl out' % sys.argv[0]
    sys.exit(1)

  vrml = Vrml()
  with open(sys.argv[1], 'r') as f:
    vrml.load(f)

  with open(sys.argv[2], 'w') as f:
    cloud = Cloud(vrml.points)
    cloud.save(f)

    

