# 16-channel custom-coil array for UIG (upper and Infra Glottic) airway MRI at 3T
## Reference paper: 
W. Alam, S. Reineke, M. R. Viswanath, R. Z. Rusho, D. V. Daele, D. Meyer, J. Liu and S. G. Lingala, "A flexible 16 channel custom coil array for accelerated imaging of upper and infra glottic airway at 3 Tesla," submitted online to Magnetic Resonance in Medicine as a Research Article.

## Paper link: 
< will be posted once online >

## Main advantages of the custom airway coil: 
- In contrast to existing head and head-neck coils, the proposed coil is in a light weight, flexible manner allowing for conformation to different face/neck anatomies, and imaging in the lateral recumbent posture.

  <img src="/images/Fig_1_caption_redo.png" alt="Alt text" width="600" title="the 16-channel custom airway coil">
- Improved SENSE parallel imaging performance for the task of static 3D imaging of the upper airway and infra glottic airway.

  <img src="/images/Fig_3_reformatted.png" alt="Alt text" width="400" title="SENSE parallel imaging performance">
- With blinded image quality ratings from two experts, we found that in the supine posture, the proposed coil enabled a net 3-fold acceleration, when accelerated either in the SI or AP direction; and a net 4 fold acceleration when accelerated in both SI and AP directions.

  <img src="/images/clubbed_along_acceleration.png" alt="Alt text" width="600" title="blinded image quality rating">
- Gain in CNR for the proposed airway coil along the entire vocal tract.

  <p float="left">
    <img src="/images/Fig_6_intsty_normalizd.png" width="400" title="CNR comparison" />
    <img src="/images/CNR_table.png" alt="Alt text", width="400", title="CNR table"/> 
  </p>
 
 ## Dataset:
  Share in-house dataset used for the in-vivo experiments once the paper is accepeted and online.
## Description of matlab and python scripts:
`SENSE_3D.m`: matlab implementation of 3D SENSE (for reconstruction)with no regularization and pseudo replica method (to evaluate g-maps)<sup>[1].

`ttGRAPPA_main.m`: matlab implementatipon of through-time spiral GRAPPA parallel imaging-based reconstruction<sup>[2].

`likert_scale_rating_violinplot.ipynb`: python scripts for demonstrating acceleration-wise violinplots for all the airway regions of interest.
## References:
- [1] Pruessmann, Klaas P., et al. "SENSE: sensitivity encoding for fast MRI." Magnetic Resonance in Medicine: An Official Journal of the International Society for Magnetic Resonance in Medicine 42.5 (1999): 952-962.
- [2] Lingala, S. G. et al. Feasibility of through-time spiral generalized autocalibrating partial parallel acquisition for low latency accelerated real-time MRI of speech. Magn. Reson. Med. 78 (6), (2017): 2275- 2282.
  
## Contact: 
  The code is meant for reproducible research. In case of any difficulty, please open an issue or directly email me at mohammadwahidul-alam@uiowa.edu
