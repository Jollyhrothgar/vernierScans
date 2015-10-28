FITTING NOTES:

Shifting Beam Width, in parallel dimension as beam displacement:

I noticed, I suppose rather obviously that, when we adjust the beam width
parameter, whose width is in the same dimension as beam displacement, that we
effectively are changing the beam offset. For example, say we have two
identically dense beams, and we scale that density profile in one dimension.
Then, we would need to displace the beams further to obtain the same general
vertex profile shape, though this shape would be spread over a wider z-vertex
range.

Therefore, for a fixed displacement, adjusting the profile which is in the same
dimension of the displacement, we can increase or decrease peak splitting. This
may be a real effect in data, and might not be easily distinguished from actual
differences between beam displacement and beam distributions.

Shifting Beam Width, in perpendicular dimension to beam displacement:

Given observation of shifting beam width in parallel dimension of beam
displacement, I expected to see no changes in the vertex-profile shape based on
the fact that the total number of ions overlapped will not change with scaling
of the vertical beam width, if the beams are being displaced horizontally. I
have confirmed this hypothesis. Scaling the beam width in the dimension
perpendicular to displacement has no effect on the vertex profile - in a brute
force computational environment, shifting the vertical beam width during the
horizontal scan will have no effect.

Shifting Beta Star
When beta-star is scaled to a smaller value, I notice the split peaks become
closer, when beta-star is scaled to a larger value, I see the right hand peak
grow, and the left hand peak shirnk. I think I am observing the combined effects
of beta-star and the crossing angle. Absent a crossing angle, the beta-star
parameter (at a value of 1) yields a narrow vertex profile centered at z=0. The
width of this distribution is narrow, about 20 cm. I am not entirely sure where
this width comes from, as the z-profile widths for the bunch in the simulation
are all larger than 20 cm. When beta-star increases, the narrow central peak
splits into two peaks, with the value of beta-star correllating to wider split
peaks. The peak width of each individual peak does not seem to shift very much
based on the value of beta-star. Generally, I can say that beta-star determines
if there is any splitting of the z-vertex profile. For large beta-star, peaks
are quite split, for small beta-star, peaks are narrowly split. The lower
beta-star the narrower the peaks seem to be, the larger, the wider.

Shifting crossing angle
Shifting the crossing angle has the affect of changing the relative peak heights
of the split peaks in the z-vertex profile. Negative crossing angle makes the
right peak higher, whereas positive crossing angle makes the left peak higher.
