# NLH: A Blind Pixel-level Non-local Method for Real-world Image Denoising
* This repository is for NLH introduced in the following paper

##### Yingkun Hou, Jun Xu, Mingxia Liu, Guanghai Liu, Li Liu, Fan Zhu, and Ling Shao, "NLH: A Blind Pixel-level Non-local Method for Real-world Image Denoising", IEEE Transactions on Image Processing, vol. 29, pp. 5121-5135, 2020.

* The code is built on Matlab 2018a.

* The folder `NLH_AWGN_Gray` includes gray scale images Additive Gaussian White Noise denosing codes and some gray scale images.

* The folder `NLH_CC` includes CC dataset and denoising codes. The related information refers to "S. Nam, Y. Hwang, Y.  Matsushita, and S. J. Kim. A holistic approach to cross-channel image noise modeling and its application to image
denoising. In CVPR, pages 1683–1691, 2016."

* The folder `NLH_DND` only includes denoising codes for DND dataset, the dataset can be downloaded from the following website:     https://noise.visinf.tu-darmstadt.de/. The related information refers to "T.  Plo¨tz  and  S.  Roth.   Bench                                 algorithms  with  real photographs.  In CVPR, 2017."
 
### Note: If one wants to conduct AWGN denoising on color images, please change "sigma_est_b = sigma_est\*4.0" to "sigma_est_b = sigma_est\*1.0" in NLH_CC.m line 47 and NLH_DND.m line 72 respectively.
