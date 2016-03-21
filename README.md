# segsalsa

This includes the source code (MATLAB) for the SegSALSA family of algorithms and some examples of its application.
The SegSALSA family of algorithms solves problems of classification with context (such as segmentation) in a convex formulation by having a hidden field driving the discrete labels.
A marginal maximum a posteriori segmentation is then performed across the discrete labels, transforming the discrete optimization problem into a continuous optimization problem.
This increases the flexbility of the prior choice.

Currently, vectorial total variation (VTV), graph total variation (GTV), and patch-based schatten norm minimization leading to structure tensor regularization (STR) are implemented.

For a more detailed overview, please refer to:

[1] J. Bioucas-Dias, F. Condessa, and J. Kovačević, "Alternating direction optimization for image segmentation using hidden Markov measure field models" (IS&T/SPIE EI 2014)

https://www.researchgate.net/publication/262986051_Alternating_direction_optimization_for_image_segmentation_using_hidden_Markov_measure_field_models

[2] F. Condessa, J. Bioucas-Dias, and J. Kovačević, "Supervised hyperspectral image segmentation: a convex formulation using hidden fields," (IEEE WHISPERS'14) 

https://www.researchgate.net/publication/299123249_SUPERVISED_HYPERSPECTRAL_IMAGE_SEGMENTATION_A_CONVEX_FORMULATION_USING_HIDDEN_FIELDS

[3] F. Condessa, J. Bioucas-Dias, and J. Kovačević, "SegSALSA-STR: a convex formulation to supervised hyperspectral image segmentation using hidden fields and structure tensor regularization," (IEEE WHISPERS'15)

https://www.researchgate.net/publication/275588151_SegSALSA-STR_A_convex_formulation_to_supervised_hyperspectral_image_segmentation_using_hidden_fields_and_structure_tensor_regularization
