# cvwrap
A Maya wrap deformer that is faster than Maya's wrap deformer, can be rebounded, has a GPU implementation, and supports inverted front of chain blend shapes.

```python
sphere = cmds.polySphere(sx=10, sy=10)[0]
cube = cmds.polyCube(w=2.3, h=2.3, d=2.3, sx=5, sy=5, sz=5)[0]

# Create a new wrap
wrap_node = cmds.cvWrap(sphere, cube, name='wrapnode', radius=0.1)

# Rebind a vertex
cmds.select(['{0}.vtx[75]'.format(sphere)])
cmds.select(['{0}.{1}'.format(cube, faces) for faces in ['f[110:111]', 'f[115:116]']], add=True)
cmds.cvWrap(rb=wrap_node)

file_path = r'E:\My Documents\maya\projects\default\data\binding.wrap'
# Export the binding
cmds.cvWrap(wrap_node, ex=file_path)

# Recreate the wrap node with the adjusted binding
cmds.delete(wrap_node)
wrap_node = cmds.cvWrap(sphere, cube, name=wrap_node, b=file_path)

# Import the binding again
cmds.cvWrap(wrap_node, im=file_path)
```
