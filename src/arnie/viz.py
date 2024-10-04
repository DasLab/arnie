###
# File is generally for functions that relate to various ways to visualize RNA sequences and experimental data.
###

import matplotlib.pyplot as plt
import matplotlib.patches as patches

def plot_structure_heatmap(data):
  # Setup constants
  dpi = 100.0
  font_size = 12.0
  padded_size = 18.0
  scaling_ratio = padded_size/dpi

  # Set default text to monospace to align sequence and structures
  plt.rcParams["font.family"] = "monospace"
  
  # Figure Height = sequence + data (experimental, predicted) + predicted structures + (opt) ROI structure
  rows = (1+len(data["reactivity"])+len(data["predictions"]))
  if "control_structure" in data.keys():
    rows += 1
  fig_height = rows * scaling_ratio

  # Figure Width determined by sequence length
  fig_width = len(data["sequence"]) * scaling_ratio

  # Create figure
  fig, axs = plt.subplots(
    nrows=rows,
    sharex=True,
    figsize=(fig_width, fig_height),
  )

  # SETUP COMMON TO ALL AXES
  # Normalize for data range from -2 to 2
  norm=plt.Normalize(-2,2)
  # Turn off the axis bounding box
  [ax.spines[:].set_visible(False) for ax in axs]
  # Set all axes x_limits from 0 to the length of the sequence
  [ax.set_xlim([0, len(data["sequence"])]) for ax in axs]
  # Hide the tick marks on the x axis for all axes
  [ax.tick_params(axis='x',length=0) for ax in axs]
  # Hide the tick marks on the y axis for all axes
  [ax.set_yticks([]) for ax in axs]
  [ax.tick_params(axis='y',length=0,pad=6) for ax in axs]
  # Plot blank image data for all axes (simplifies alignment of text in sequence and prediction plotting)
  [ax.imshow([[0]*(len(data["sequence"])+1)], cmap="bwr",norm=norm) for ax in axs]
  # Collapse sub-plot spacing
  fig.subplots_adjust(hspace=0)
  
  # Set the plot title
  axs[0].set_title(data["title"], fontweight='bold')

  ################################################
  # PLOTTING
  ################################################
  ax_index = 0
  
  # Plot sequence text
  axs[0].set_yticks([0],labels=["sequence"])
  for (i, char) in enumerate(data["sequence"]):    
    if (i>=0):
      axs[0].text(i,0,char,fontfamily='monospace', ha="center", va="center")

  # Create the reactivity data heatmap
  for (j, data_label) in enumerate(data["reactivity"]):
    ax_index += 1 
    axs[ax_index].set_yticks([0],labels=[data_label])
    reactivity = data["reactivity"][data_label]["data"]
    BLANK_OUT5 = data["reactivity"][data_label]["BLANK_OUT5"]
    BLANK_OUT3 = data["reactivity"][data_label]["BLANK_OUT3"]

    # Reactivity needs to be a list of numbers
    if type(reactivity) != list or type(reactivity[0]) != float:
      print("WARNING: reactivity data in unexpected format")
        
    # If the data has blank out regions, add them.
    display_data = [-1.0] * BLANK_OUT5 + reactivity + [-1.0] * BLANK_OUT3

    # Plot the heatmap and the blank out regions
    pos = axs[ax_index].imshow([display_data], cmap="bwr", norm=norm)
    blank5 = patches.Rectangle((-0.5,-0.5),BLANK_OUT5,1,color="gray")
    axs[ax_index].add_patch(blank5)
    blank3 = patches.Rectangle((len(display_data)-BLANK_OUT3-0.5,-0.5),BLANK_OUT3,1,color="gray")
    axs[ax_index].add_patch(blank3)
      
  plt.colorbar(pos, ax=axs)

  # Plot control structure (if provided)
  if "control_structure" in data.keys():
    if type(data["control_structure"]["start_index"]) != float and type(data["control_structure"]["structure"]) != str:
      return
    ax_index += 1
    axs[ax_index].set_yticks([0],labels=["Control Structure"])
    start_index = data["control_structure"]["start_index"]
    for (i, char) in enumerate(data["control_structure"]["structure"]):    
      if (i>=0):
          axs[ax_index].text(i+start_index,0,char,fontfamily='monospace', ha="center", va="center")
            
  # Plot predictions text
  for (j, predictor_name) in enumerate(data["predictions"]):
    ax_index += 1
    axs[ax_index].set_yticks([0],labels=[predictor_name])
    for (i, char) in enumerate(data["predictions"][predictor_name]):    
      if (i>=0):
        axs[ax_index].text(i,0,char,fontfamily='monospace', ha="center", va="center")
    if ax_index == rows-1:
      axs[ax_index].tick_params(axis='x',length=4,direction='out')
          
  return fig