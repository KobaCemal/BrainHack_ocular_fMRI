# libraries
from ipywidgets import interact, IntSlider, FloatSlider
from matplotlib.colors import LinearSegmentedColormap


def plot_mri(image_data, projection: str, mask_data=None, *, initial_channel=0, initial_time=0, mask_color='Wistia'):
    '''
    image_data: takes 4d matrix as input\n\n
    projection: takes XY, XZ, YZ as arguments\n\n
    mask_data: takes 4d matrix as input \n
    IMPORTANT: same size as image_data\n
    mask_color: takes cmap as argument 
    '''


    # adjusting output
    match projection:
        case 'XY':
            channel_max = image_data.shape[2]
            channel_name = 'Z'
        case 'XZ':
            channel_max = image_data.shape[1]
            channel_name = 'Y'
        case 'YZ':
            channel_max = image_data.shape[0]
            channel_name = 'X'

    # making mask transparent 
    if mask_data is not None:        
        # making zeros transparent 
        mask_data = np.ma.masked_where(mask_data == 0, mask_data)
        # printing shapes of images (to make sure nothing is wrong)
        print(f'Image shape: {image_data.shape}')
        print(f'Mask shape:  {mask_data.shape}')

    # plot part 
    def explore_img(Time, channel, alpha):
        # fig settings
        plt.figure(figsize=(12, 8))
        # printing 
        match projection:
            case 'XY':
                plt.imshow(image_data[:, :, channel, Time], cmap='gray')    
                if mask_data is not None:        
                    plt.imshow(mask_data[:, :, channel, Time], cmap=mask_color, alpha=alpha)       
            case 'XZ':
                plt.imshow(image_data[:, channel, :, Time], cmap='gray')
                if mask_data is not None:        
                    plt.imshow(mask_data[:, channel, :, Time], cmap=mask_color, alpha=alpha)  
            case 'YZ':
                plt.imshow(image_data[channel, :, :, Time], cmap='gray')
                if mask_data is not None:        
                    plt.imshow(mask_data[channel, :, :, Time], cmap=mask_color, alpha=alpha)  
        
        plt.axis('off')
        # plt.show()

    # interactive function
    interact(explore_img, 
        Time=IntSlider(min=0, max=image_data.shape[3]-1, step=1, value=initial_time),
        channel=IntSlider(min=0, max=channel_max -1, step=1, value=initial_channel, description=channel_name),
        alpha=FloatSlider(min=0, max=1, step=0.01, value=0.2, description='Transparency')
        )

