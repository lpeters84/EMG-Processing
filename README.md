# EMG-Processing
Collection and MATLAB processing of electromyography data  
   
![EMG_signals](https://github.com/user-attachments/assets/255faf6d-03b7-45bb-aa53-d2778dfc8813)


# Description
This data was collected using a BioRadio and surface EMG electrodes on the biceps brachii. 5 contractions with rest in between were performed, each aiming to increase the magnitude of each subsequent contraction. 
  
EMG_filtered was prefiltered by the Bioradio.  
EMG_raw is the raw data collected by the BioRadio. 
  
This code processes the raw data and compares it to the prefiltered data.

The timing of each muscle contraction and the corresponding integrated contraction values are then calculated to estimate EMG values.
