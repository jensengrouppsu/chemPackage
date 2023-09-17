from __future__ import division, print_function

def fit_func(function, parameters, xdata, ydata, samples=1000):
    '''This fits data using SciPy's leastsq routine, but then uses the
    bootstrap method to find the fit quality.  This essentially performes a
    monte carlo analysis on the function using the fitted parameters, and
    determines the more correct parameters and confidence based on that
    analysis.

    function should be a function where the first argument is the x-values,
    and the second argument is a list of the constants of the function (that
    you wish to fit).

    parameters is a list of the initial guesses of the constants you wish to
    fit.  The order must be the same order accepted in the function.

    xdata is the raw x values

    ydata is the raw y values

    samples is the number of monte carlo sample points to take.
    The default is 1000.
    '''
    from scipy.optimize import leastsq
    from numpy import array, mean, std
    from numpy.random import normal

    # First, define a factory function for the error of the data
    def error(params, x, y):
        return y - function(x, params)

    # Define a fitting function
    def fit_helper(e, p, x, y):
        ret = leastsq(e, p, args=(x, y), full_output=True)
        if ret[4] not in set([1, 2, 3, 4]):
            raise ValueError (ret[3])
        else:
            return ret[0]

    # Now use the above function to send data to least squares
    params = fit_helper(error, parameters, xdata, ydata)

    # Get the residuals of this fit
    res   = function(xdata, params) - ydata
    s_res = std(res)

    # Use a monte carlo technique to find error approximations
    # We basically are varying the y data with errors sampled on a normal
    # distribution, and finding the fits with that new data
    new_params = []
    for i in xrange(samples):
        new_y = ydata + normal(0, s_res, len(ydata))
        new_params.append(fit_helper(error, parameters, xdata, new_y))
    new_params = array(new_params)

    # Now find the mean and st. dev. of the new parameters
    return mean(new_params, 0), std(new_params, 0)

##########################################################################
# NOTE: The Parameter class and the old_fit function go together, but they
# should be considered old and therefore should be avoided.
##########################################################################

class Parameter:
    '''This class defines the parameters (read: variables) used for
    fitting.  The parameters should be set with an initial guess.
    '''
    def __init__(self, value):
        self.value = value

    def set(self, value):
        self.value = value

    def __call__(self):
        return self.value

def fit(function, parameters, y, x = None):
    '''This function receives the function to fit (using the parameters
    defined above), the parameters, and the raw data (where the domain is
    assumed to be the whole integers corresponding to your y values if
    x is omitted).  Returns the parameters with the correct values.
    ''' 
    from scipy import optimize
    from numpy import arange
    def f(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        return y - function(x)

    if x is None: x = arange(y.shape[0])
    p = [param() for param in parameters]
    params, sucess = optimize.leastsq(f, p)
    if sucess not in set([1, 2, 3, 4]): raise ValueError ('No fit found')
    return params

def tick_scale(minval, maxval, numticks=None, reverse=False):
    '''This function will find the best points for tick marks on a scale
    bar given some range of data.  It will shoot for as close to seven
    ticks as possible without going under.  This can be overriden by
    choosing a specific number of ticks.  A list of the tick locations
    is returned.

    KNOWN BUG: This function only works for large values.
    '''
    from prep import range_check
    from mfunc import roundint

    # Make sure low is lower than high
    low, high = range_check(minval, maxval)

    # Find the difference in the data
    diff = maxval - minval

    # Try different multiples of 5 until we find one that gives close
    # to 9 ticks.  Only do this if numticks was not defined.
    if numticks is None:
        multiple = 5
        while True:
            ticks = int(diff / multiple)
            if ticks <= 9:
                break
            else:
                multiple += 5
    # Round the skipping interval to the nearest multiple of five if numticks
    # is given.
    else:
        multiple = intround(diff / numticks, 5)

    # Find appropriate endpoints
    lowend = intround(low, multiple)
    if lowend < int(low): lowend =+ multiple
    highend = intround(high, multiple)
    if highend > int(high): highend = highend - multiple
    # Return the range. 
    if reverse:
        return reversed(range(lowend, highend+multiple, multiple))
    else:
        return range(lowend, highend+multiple, multiple)
