# InSAR-phase-jump-correction
An algorithm to automatically detect and correct phase jumps in unwrapped InSAR time series

This code comes along with the paper: Fanghui Deng, Mark Zumberge. Seafloor motion from offshore platforms using satellite radar images – a case study in the Adriatic Sea. Submitted. 

The code was written and tested in Python 3.8.10 under Ubuntu 20.04. The following are steps of how to use the code based on the example time series ts_example.txt

1. Run step1_rotate_time_series_with_plots.py. This will rotate the time series using different angles, and calculate corresponding histograms and density curves. If you do not want the output plots, you can run step1_rotate_time_series_without_plots.py instead. This will speed up the process. 
2. Run step2_decide_best_rotation.py. This will find the best rotation angle. 
3. Run step3_correct_ts_basedon_best_rotation.py. This will correct the original time series based on the best rotation angle.

For questions and feedbacks, please contact Fanghui Deng at fadeng@ucsd.edu.

This is research code provided to you "AS IS" with no warranties of correctness. Use at your own risk.
