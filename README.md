# hemidentification

This is supporting code for a project designed to identify whether a
partial within-hemisphere connectome belongs to a left or a right hemisphere
using Human Connectome Project-Young Adult (HCP-YA) data
(Van Essen et al., 2013).

The project especially relies on results and code from Hannum et al. (2023).

## Description of code

`code/fake_hemiconnectome/`: Code to generate a schematic (hemi-)connectome
    with four symmetric ROIs per hemisphere, a la Glasser parcellation.

`code/fold_testing/`: Code to test the creation of *k*-folds for model
    evaluation, using HCP-YA demographic data.

`code/formula`: TeX code to generate the formula for Matthew's correlation
    coefficient, because the Google Docs equation editor is inadequate.

`code/hypotheses/`: Code to generate the hypothesis figures from the registered
    report.

## References

 - Glasser, M. F., Coalson, T. S., Robinson, E. C., Hacker, C. D., Harwell, J.,
    Yacoub, E., Ugurbil, K., Andersson, J., Beckmann, C. F., Jenkinson, M.,
    Smith, S. M., & Van Essen, D. C. (2016).
    A multi-modal parcellation of human cerebral cortex. *Nature*, 536(7615),
    Article 7615. https://doi.org/10.1038/nature18933

 - Hannum, A., Lopez, M. A., Blanco, S. A., & Betzel, R. F. (2023).
    High-accuracy machine learning techniques for functional connectome
    fingerprinting and cognitive state decoding. *Human Brain Mapping*, 44(16),
    5294–5308. https://doi.org/10.1002/hbm.26423

 - Van Essen, D. C., Smith, S. M., Barch, D. M., Behrens, T. E. J., Yacoub, E.,
    & Ugurbil, K. (2013). The WU-Minn Human Connectome Project: An overview.
    *NeuroImage*, 80, 62–79. https://doi.org/10.1016/j.neuroimage.2013.05.041

