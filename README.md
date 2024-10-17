# Misc_code_public

Misc. code that doesn't fall into other repositories, to be publicly accessible.

### count_cells_stereoscope_Ngo2024.m (MATLAB)

Read 2D stereoscope images, detect cells by thresholding,  culling small and large objects. (Objects = immune cells.)
Then high-pass filter and again threshold, etc.; combine (OR) the two binary images.

For data and analysis described in Ngo et al. 2024 . Copied "as used" with code dependent on particular file structures. However: the function findObjects() takes an image as an input, and is generally applicable, regardless of file structure.

More importantly, the code will illustrate what analysis operations were performed.

