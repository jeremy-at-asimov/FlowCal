import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import FlowCal.gate

default_stats = ['mean', 'gmean', 'mode', 'median', 'std', 'cv', 'gstd', 'gcv', 'iqr', 'rcv']
default_bin_edges = [-100, -50, 0, 50, 1e2, 3.16e2, 1e3, 3.16e3, 1e4, 3.16e4, 1e5, 3.16e5, 1e6, 3.16e6, 1e7]
default_diag_edges = [-100, 1e2, 3.16e2, 1e3, 3.16e3, 1e4, 3.16e4, 1e5, 3.16e5, 1e6, 3.16e6, 1e7]

def bin1d(bin_by, list_FCSData, bin_edges=default_bin_edges, stats=default_stats, channels=[], events_threshold=0):
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
    events_threshold : int, optional
        events needed for stats to be calculated for a given bin

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

                    if len(within_bin) < events_threshold:
                        temp_df.loc[:] = np.nan

                    row_df = pd.concat([row_df, temp_df], axis=1, sort=False, join='outer')

                sample_df = pd.concat([sample_df, row_df], axis=0, sort=False, join='outer')

        df = pd.concat([df, sample_df], axis=1, sort=False, join='outer')

    return (df)

def bin_diag(list_FCSData, bin_edges=default_diag_edges, stats=default_stats, channels=[0,1], events_threshold=0, xscale='logicle', yscale='logicle'):
    """
    Bins a list of FCSData objects and generates a dataframe with
    statistics calculated within each bin. Bins are diagonally
    slanted trapezoids rather than rectangular, which may better
    account for different samples with high signal different axes.

    Parameters
    ----------
    list_FCSData : a list of FCSData objects
        NxD flow cytometry data where N is the number of events and D is
        the number of parameters (aka channels).
    bin_diag_edges : list of int
        Values (often fluorescence) within channel used to define the
        trapezoidal bins.
    stats : list of str
        names for statistics to be calculated. Values include ['mean',
        'gmean', 'mode', 'median', 'std', 'cv', 'gstd', 'gcv', 'iqr',
        'rcv'].
    channels : list of str or int
        list of channels to make statistics for
    events_threshold : int, optional
        events needed for stats to be calculated for a given bin
    xscale : str, optional
        Scale for x-axis when determining if datapoint are within a bin
    yscale : str, optional
        Scale for y-axis when determining if datapoint are within a bin

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

        sample_df = pd.DataFrame()
        contours = []

        min_edge = min(bin_edges)
        max_edge = max(bin_edges)

        # make trapezoidal gating contours for low expression region
        for i, (low_edge, high_edge) in enumerate(zip(bin_edges[:-1], bin_edges[1:])):
            contour = [(high_edge,min_edge),
                       (min_edge,high_edge),
                       (min_edge,low_edge),
                       (low_edge,min_edge)]
            if contour[-1] != contour[0]: # close the contour if necessary
                contour.append(contour[0])
            contour = np.asarray(contour)
            contours.append(contour)

        # make trapezoidal gating contours for high expression region
        for i, (high_edge, low_edge) in enumerate(zip(reversed(bin_edges[:-1]), reversed(bin_edges[1:]))):
            contour = [(max_edge, high_edge),
                       (high_edge, max_edge),
                       (max_edge, low_edge),
                       (low_edge, max_edge)]
            if contour[-1] != contour[0]:  # close the contour if necessary
                contour.append(contour[0])
            contour = np.asarray(contour)
            contours.append(contour)

        # loop over trapezoidal bins and calculate statistics
        for i, contour in enumerate(contours):
            contour = [contour]
            within_bin = FlowCal.gate.contour2d(temp_FCSData, channels=channels, gate_contour=contour,
                                               xscale=xscale, yscale=yscale, full_output=False)

            if within_bin.shape[0] > 0:
                row_df = pd.DataFrame()

                for stat in stats:
                    function = getattr(FlowCal.stats, stat)
                    temp_np_array = np.copy(empty_np)
                    # print(within_bin.shape)
                    temp_np_array[:] = function(within_bin)

                    temp_df = pd.DataFrame(data=temp_np_array, index=['bin_'+str(i)], columns=channel_names)
                    temp_df = temp_df[channels]
                    filtered_channel_names = list(temp_df.columns)
                    temp_df.columns = [str(temp_FCSData) + '_' + stat + '_' + channel_name for channel_name in
                                       filtered_channel_names]

                    if len(within_bin) < events_threshold:
                        temp_df.loc[:] = np.nan

                    row_df = pd.concat([row_df, temp_df], axis=1, sort=False, join='outer')

                sample_df = pd.concat([sample_df, row_df], axis=0, sort=False, join='outer')

        df = pd.concat([df, sample_df], axis=1, sort=False, join='outer')

    return (df)

class model_wrapper:

    def __init__(self):
        pass

    def linear_model(self, x, m, b):
        if self.transform == 'log':
            return np.log(m * np.exp(x) + b)
        elif self.transform == 'logicle':
            x_logicle = FlowCal.plot._LogicleTransform(data=x, channel=0)
            x_inverse = FlowCal.plot._InterpolatedInverseTransform(transform=x_logicle, smin=0, smax=x_logicle.M)
            return x_inverse.transform_non_affine(m * x_logicle.transform_non_affine(x) + b)
        else:
            return m * x + b

'''
def fit_linear(x, m, b):
    return m * x + b


def fit_log(x, m, b):
    return np.log(m * np.exp(x) + b)


def fit_logicle(x, m, b):
    x_logicle = FlowCal.plot._LogicleTransform(data=x, channel=0)
    x_inverse = FlowCal.plot._InterpolatedInverseTransform(transform=x_logicle, smin=0, smax=x_logicle.M)

    return x_inverse.transform_non_affine(m * x_logicle.transform_non_affine(x) + b)
'''

def model_fit(df,
              xaxis='',
              yaxis='',
              transform='log',
              model='default'):
    """
    Generate model fits from binned data in pandas dataframe.
    Added by Jeremy Gam, jeremy@asimov.io

    Parameters
    ----------
    df : pandas dataframe
        Binnedlow cytometry data to plot.
    xaxis : str
        Name of statistic and channel for fit on x axis (e.g.
        'median_PB450-A').
    yaxis : str
        Name of statistic and channel for fit on y axis (e.g.
        'median_APC-A').

    Notes
    -----
    Fits with logicle transform seem sensitive to noise so bins
    with low numbers of events should be excluded in steps prior
    to fitting

    """
    # Check appropriate inputs
    if not xaxis:
        raise ValueError('need string for x axis statistic + channel')
    if not yaxis:
        raise ValueError('need string for y axis statistic + channel')

    # Get sample names from df
    sample_names = []
    for column in df.columns:
        xaxis_index = column.find('_' + xaxis)
        if xaxis_index != -1:
            sample_names.append(column[0:xaxis_index])

    # Iterate through sample names. Transform and fit data
    fits = []
    for i, sample_name in enumerate(sample_names):
        xcolumn = sample_name + '_' + xaxis
        ycolumn = sample_name + '_' + yaxis
        xdata = df[xcolumn].values
        ydata = df[ycolumn].values

        to_keep = ~np.isnan(xdata) & ~np.isnan(ydata)  # remove nans prior to fits
        xdata = xdata[to_keep]
        ydata = ydata[to_keep]

        if transform == 'logicle':

            x_logicle = FlowCal.plot._LogicleTransform(data=xdata, channel=0)
            x_inverse = FlowCal.plot._InterpolatedInverseTransform(transform=x_logicle, smin=0, smax=7)  # x_logicle.M)

            xdata_transform = x_inverse.transform_non_affine(xdata)
            ydata_transform = x_inverse.transform_non_affine(ydata)

        elif transform == 'log':

            xdata_log = np.log(xdata)
            ydata_log = np.log(ydata)

            to_keep = ~np.isnan(xdata_log) & ~np.isnan(ydata_log)  # remove zeros prior to log transforms

            xdata_transform = xdata_log[to_keep]
            ydata_transform = ydata_log[to_keep]

        else:

            xdata_transform = xdata
            ydata_transform = ydata


        if model == 'hill':
            pass #for now

        elif model == 'linear':
            model_class = model_wrapper()
            model_class.transform = transform
            popt, pcov = curve_fit(model_class.linear_model, xdata_transform, ydata_transform)

        elif model == 'default':
            print('Warning: no model given, using linear model as default')
            model_class = model_wrapper()
            model_class.transform = transform
            popt, pcov = curve_fit(model_class.linear_model, xdata_transform, ydata_transform)

        else:
            print('Error: given model not implemented yet, no fit will be provided')
            popt = 0

        fits.append(popt)

    return fits