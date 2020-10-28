import numpy as np
import pandas as pd
import FlowCal.gate

default_stats = ['mean', 'gmean', 'mode', 'median', 'std', 'cv', 'gstd', 'gcv', 'iqr', 'rcv']
default_bin_edges = [-100, -50, 0, 50, 1e2, 3.16e2, 1e3, 3.16e3, 1e4, 3.16e4, 1e5, 3.16e5, 1e6, 3.16e6, 1e7]

def bin1d(bin_by, list_FCSData, bin_edges=default_bin_edges, stats=default_stats, channels=[]):
    """
    Bins a list of FCSData objects and generates a dataframe with
    statistics calculated within each bin.

    Parameters
    ----------
    bin_by : str
        channel name used to generate bins.
    list_FCSData : a list of FCSData objects
        NxD flow cytometry data where N is the number of events and D is
        the number of parameters (aka channels).
    bin_edges : list of int
        Values (often fluorescence) within channel used to divide bins.
    stats : list of str
        names for statistics to be calculated. Values include ['mean',
        'gmean', 'mode', 'median', 'std', 'cv', 'gstd', 'gcv', 'iqr',
        'rcv'].
    channels : list of str
        list of channels to make statistics for

    Returns
    -------
    pandas dataframe
        df with column names including the sample name, statistic, and
        channel.

    Notes
    -----


    """
    first_FCSData = list_FCSData[0]
    empty_np = np.empty((1, first_FCSData.shape[1]))
    empty_np[:] = np.NaN
    first_FCSData = list_FCSData[0]
    channel_names = list(first_FCSData.channels)
    df = pd.DataFrame()

    if not channels:
        channels = channel_names

    # iterate over all FCSData objects
    for temp_FCSData in list_FCSData:

        # iterate over bin edges and find stats
        sample_df = pd.DataFrame()
        for i, (left_edge, right_edge) in enumerate(zip(bin_edges[:-1], bin_edges[1:])):

            within_bin = FlowCal.gate.high_low(temp_FCSData, channels=bin_by, low=left_edge, high=right_edge,
                                               full_output=False)

            if within_bin.shape[0] > 0:
                row_df = pd.DataFrame()

                for stat in stats:
                    function = getattr(FlowCal.stats, stat)
                    temp_np_array = np.copy(empty_np)
                    # print(within_bin.shape)
                    temp_np_array[:] = function(within_bin)
                    temp_df = pd.DataFrame(data=temp_np_array, index=[str(left_edge)], columns=channel_names)
                    temp_df = temp_df[channels]
                    filtered_channel_names = list(temp_df.columns)
                    temp_df.columns = [str(temp_FCSData) + '_' + stat + '_' + channel_name for channel_name in
                                       filtered_channel_names]

                    row_df = pd.concat([row_df, temp_df], axis=1, sort=False, join='outer')

                sample_df = pd.concat([sample_df, row_df], axis=0, sort=False, join='outer')

        df = pd.concat([df, sample_df], axis=1, sort=False, join='outer')

    return (df)