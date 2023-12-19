"""Utilities for the GRSS orbit propagation code"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

__all__ = [ 'equat2eclip',
            'eclip2equat',
            'plot_solar_system',
            'plot_ca_summary',
            'plot_bplane'
]

earth_obliq = 84381.448/3600.0*np.pi/180.0
equat2eclip = np.array([[1.0,          0.0,                 0.0        ],
                        [0.0,  np.cos(earth_obliq), np.sin(earth_obliq)],
                        [0.0, -np.sin(earth_obliq), np.cos(earth_obliq)]])
eclip2equat = np.array([[1.0,          0.0,                 0.0        ],
                        [0.0,  np.cos(earth_obliq), -np.sin(earth_obliq)],
                        [0.0,  np.sin(earth_obliq),  np.cos(earth_obliq)]])

def _add_planet_xy(axis, dist, color, alpha, lstyle):
    """
    Add a planet to the x-y plane plot.

    Parameters
    ----------
    axis : matplotlib axis
        axis to plot on
    dist : float
        semi-major axis of the planet's orbit
    color : str
        color of the planet's orbit
    alpha : float
        transparency of the planet's orbit
    lstyle : str
        linestyle of the planet's orbit

    Returns
    -------
    None : NoneType
        None
    """
    circle = plt.Circle((0, 0), dist, color=color, fill=False, alpha=alpha, linestyle=lstyle)
    axis.add_artist(circle)
    return None

def _add_planet_xz(axis, dist, obliq, color, alpha, lstyle):
    """
    Add a planet to the x-z plane plot.

    Parameters
    ----------
    axis : matplotlib axis
        axis to plot on
    dist : float
        semi-major axis of the planet's orbit
    obliq : float
        obliquity of the planet's orbit
    color : str
        color of the planet's orbit
    alpha : float
        transparency of the planet's orbit
    lstyle : str
        linestyle of the planet's orbit

    Returns
    -------
    None : NoneType
        None
    """
    sin = dist*np.sin(obliq*np.pi/180)
    cos = dist*np.cos(obliq*np.pi/180)
    axis.plot([-cos, cos], [-sin, sin], color=color, alpha=alpha, linestyle=lstyle)
    return None

def plot_solar_system(axis, xy_plane, alpha=0.25, lstyle='--'):
    """
    Plots the solar system in the x-y or x-z plane. For visual purposes only.

    Parameters
    ----------
    axis : matplotlib axis
        axis to plot on
    xy_plane : bool
        True to plot in x-y plane, False to plot in x-z plane
    alpha : float, optional
        transparency of the planet's orbit, by default 0.25
    lstyle : str, optional
        linestyle of the planet's orbit, by default '--'

    Returns
    -------
    None : NoneType
        None
    """
    # draw the sun
    axis.scatter(0, 0, s=50, c='yellow', marker='o', edgecolors='gray', alpha=alpha)
    # draw the orbit of the planets as circles for x-y plane
    num_planets = 9
    dists = [0.387, 0.723, 1.0, 1.524, 5.203, 9.539, 19.18, 30.06, 39.53]
    colors = ['gray', 'orange', 'blue', 'red', 'brown', 'orange', 'cyan', 'blue', 'gray']
    if xy_plane:
        for i in range(num_planets):
            _add_planet_xy(axis, dists[i], colors[i], alpha, lstyle)
    # draw the orbit of the planets as inclined lines for x-z plane
    else:
        obliqs = [7.25, 3.39, 0.0, 1.85, 1.3, 2.49, 0.77, 1.77, 0.0]
        for i in range(num_planets):
            _add_planet_xz(axis, dists[i], obliqs[i], colors[i], alpha, lstyle)
    return None

def get_scale_factor(body_id):
    """
    Get the scale factor for plotting close approaches.

    Parameters
    ----------
    body_id : int/str
        SPICE ID of the body or name of the body
    
    Returns
    -------
    scale_factor : float
        scale factor for plotting close approaches, in km
    units : str
        units of the scale factor
    """
    if body_id in {1,199,"Mercury Barycenter"}:
        scale_factor = 2440.53
        units = "R$_{Mercury}$"
    elif body_id in {2,299,"Venus Barycenter"}:
        scale_factor = 6051.8
        units = "R$_{Venus}$"
    elif body_id in {399,"Earth"}:
        scale_factor = 6378.137
        units = r"R$_\oplus$"
    elif body_id in {4,499,"Mars Barycenter"}:
        scale_factor = 3396.19
        units = "R$_{Mars}$"
    else:
        raise ValueError("Unknown body ID")
    return scale_factor, units

def plot_ca_summary(prop_sim, flyby_body, central_body='Earth',
                    x_min=None, x_max=None, y_min=1e-1, y_max=None):
    # sourcery skip: low-code-quality
    """
    Plot the relative radial velocity and distance of the close approaches of
    a given prop_sim object.

    Parameters
    ----------
    prop_sim : PropSimulation
        PropSimulation object to show the close approaches of
    flyby_body : str
        Name of the body to plot the close approaches of
    central_body : str, optional
        Name of the body to plot the close approaches with, by default 'Earth'.
        (should match name in prop_sim)
    x_min : float, optional
        Minimum time to show, by default None
    x_max : float, optional
        Maximum time to show, by default None
    y_min : float, optional
        Minimum distance to show, by default 1e-1 body radii
    y_max : float, optional
        Maximum distance to show, by default None

    Returns
    -------
    None : NoneType
        None
    """
    rel_state = np.zeros((len(prop_sim.interpParams.tStack), 6))
    t_stack = np.array(prop_sim.interpParams.tStack)
    x_stack = np.array(prop_sim.interpParams.xIntegStack)
    for idx, time in enumerate(t_stack):
        central_body_state = prop_sim.get_spiceBody_state(time, central_body)
        ca_body_state = x_stack[idx]
        rel_state[idx,0] = ca_body_state[0] - central_body_state[0]
        rel_state[idx,1] = ca_body_state[1] - central_body_state[1]
        rel_state[idx,2] = ca_body_state[2] - central_body_state[2]
        rel_state[idx,3] = ca_body_state[3] - central_body_state[3]
        rel_state[idx,4] = ca_body_state[4] - central_body_state[4]
        rel_state[idx,5] = ca_body_state[5] - central_body_state[5]
    rel_dist = np.linalg.norm(rel_state[:,:3], axis=1)
    rel_radial_vel = np.sum(rel_state[:,:3]*rel_state[:,3:6], axis=1) / rel_dist
    scale_factor, units = get_scale_factor(central_body)
    fac = prop_sim.consts.du2m/1000/scale_factor

    ca_times = [
        approach.t
        for approach in prop_sim.caParams
        if approach.flybyBody == flyby_body
        and approach.centralBody == central_body
        and approach.impact is False
    ]
    impact_times = [
        approach.t
        for approach in prop_sim.impactParams
        if approach.flybyBody == flyby_body
        and approach.centralBody == central_body
        and approach.impact is True
    ]
    ca_times = Time(ca_times, format='mjd', scale='tdb').tdb.datetime
    impact_times = Time(impact_times, format='mjd', scale='tdb').tdb.datetime
    if len(ca_times) == 0 and len(impact_times) == 0:
        print("No close approaches or impacts found")
        return None
    times = Time(t_stack, format='mjd', scale='tdb').tdb.datetime
    lwidth = 1
    fig = plt.figure(figsize=(6, 6), dpi=150)
    ax1 = plt.gca()
    ax1.plot(times, rel_radial_vel, '-', lw=lwidth, color='C0', label="Relative Radial Velocity")
    ax1.set_xlabel("Time [TDB]")
    ax1.set_ylabel("Relative Radial Velocity [au/day]", color='C0')
    for time in ca_times:
        ax1.axvline(time, ls=':', color='gray', lw=0.75)
    if len(ca_times) > 0:
        ax1.axvline(np.inf, ls=':', color='gray', lw=0.75, label="CA Time")
    for time in impact_times:
        ax1.axvline(time, ls=':', color='red', lw=0.75)
    if len(impact_times) > 0:
        ax1.axvline(np.inf, ls=':', color='red', lw=0.75, label="Impact Time")
    if x_min is not None:
        ax1.set_xlim(xmin=x_min.tdb.datetime)
    if x_max is not None:
        ax1.set_xlim(xmax=x_max.tdb.datetime)
    ax1.tick_params(axis='y', colors='C0', which='both', direction='in')
    ax2 = plt.gca().twinx()
    ax2.semilogy(times, fac*rel_dist, '--', lw=lwidth, color='C1', label="Relative Distance")
    ax2.axhline(1, ls=':', color='C1', lw=0.75, label=f"1 {units}")
    ax2.set_ylabel(fr"Relative Distance [{units}]", color='C1')
    if y_min is not None:
        ax2.set_ylim(ymin=y_min)
    if y_max is not None:
        ax2.set_ylim(ymax=y_max)
    ax2.tick_params(axis='y', colors='C1', which='both', direction='in')
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')
    fig.autofmt_xdate()
    plt.show()
    return None

def data_to_ellipse(x_data, y_data, n_std, plot_offset, bplane_type,
                    print_ellipse_params, units, sigma_points):
    """
    Convert two sets of points to an ellipse.

    Parameters
    ----------
    x_data : np.ndarray
        x-coordinates of the points
    y_data : np.ndarray
        y-coordinates of the points
    n_std : float, optional
        Number of standard deviations to plot
    plot_offset : bool, optional
        True to plot the offset from the mean, False to plot the raw data
    bplane_type : str
        type of B-plane
    print_ellipse_params : bool, optional
        True to print the ellipse parameters, False to not print
    units : str
        units of the plot
    sigma_points : prop.SigmaPoints, optional
        SigmaPoints object for reconstructing mean and cov, by default None.

    Returns
    -------
    ellipse : np.ndarray
        ellipse points
    """
    if sigma_points is None:
        x_mean = np.mean(x_data)
        y_mean = np.mean(y_data)
        cov = np.cov(x_data, y_data)
    else:
        xy_data = np.vstack((x_data, y_data)).T
        mean, cov = sigma_points.reconstruct(xy_data)
        x_mean = mean[0]
        y_mean = mean[1]
    if plot_offset:
        x_data -= x_mean
        y_data -= y_mean
    eigvals, eigvecs = np.linalg.eig(cov)
    idx = eigvals.argsort()[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:,idx]
    theta = np.arctan(eigvecs[1,0]/eigvecs[0,0])
    if theta < 0:
        theta += np.pi
    sma = np.sqrt(eigvals[0]) * n_std
    smi = np.sqrt(eigvals[1]) * n_std
    if print_ellipse_params:
        print(f'{bplane_type} ellipse sma: {sma/n_std} {units}')
        print(f'{bplane_type} ellipse smi: {smi/n_std} {units}')
        print(f'{bplane_type} ellipse theta: {theta*180/np.pi} deg')
    theta_arr = np.linspace(0, 2*np.pi, 100)
    ellipse = np.array([sma*np.cos(theta_arr), smi*np.sin(theta_arr)])
    rot_mat = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    ellipse = rot_mat @ ellipse
    if not plot_offset:
        ellipse[0,:] += x_mean
        ellipse[1,:] += y_mean
    return ellipse

def days_to_dhms(days):
    """
    Convert days to days, hours, minutes, seconds, and a summary string.

    Parameters
    ----------
    days : float
        Decimal days to convert

    Returns
    -------
    day : int
        Days
    hour : int
        Hours
    mins : int
        Minutes
    sec : float
        Seconds
    string : str
        Summary string
    """
    days = np.abs(days)
    day = int(days)
    hour = int((days - day)*24)
    mins = int(((days - day)*24 - hour)*60)
    sec = (((days - day)*24 - hour)*60 - mins)*60
    string = f'{day}d {hour:02d}:{mins:02d}:{sec:06.3f}'
    return day, hour, mins, sec, string

def plot_bplane(ca_list, plot_offset=False, scale_coords=False, n_std=3, units_km=False,
                equal_axis=True, print_ellipse_params=False, show_central_body=True,
                sigma_points=None):
    """
    Plot the B-planes of a list of close approaches.

    Parameters
    ----------
    ca_list : list
        List of CloseApproach objects
    plot_offset : bool, optional
        True to plot the offset from the mean, False to plot the raw data, by default False.
    scale_coords : bool, optional
        True to scale the coordinates to the body's radius, False to leave in km, by default False.
    n_std : float, optional
        Number of standard deviations to plot, by default 3.0
    units_km : bool, optional
        True to plot in km, False to plot in AU, by default False.
    equal_axis : bool, optional
        True to make the x and y axes equal, False to leave them as is, by default True.
    print_ellipse_params : bool, optional
        True to print the ellipse parameters, False to not print, by default False.
    show_central_body : bool, optional
        True to show the central body, False to not show it, by default True.
    sigma_points : prop.SigmaPoints, optional
        SigmaPoints object for reconstructing mean and cov, by default None.

    Returns
    -------
    None : NoneType
        None
    """
    if len(ca_list) == 0:
        print("No close approaches given (list is empty)")
        return None
    au2units = 149597870.7 if units_km else 1.0
    body_id = ca_list[0].centralBodySpiceId
    central_body_radius, units = get_scale_factor(body_id)
    if not scale_coords:
        units = "km"
    if not units_km:
        central_body_radius /= 149597870.7
        units = "AU"
    times = np.array([approach.t for approach in ca_list])
    map_times = np.array([approach.tMap for approach in ca_list])
    kizner_x = np.array([approach.kizner.x*au2units for approach in ca_list])
    kizner_y = np.array([approach.kizner.y*au2units for approach in ca_list])
    opik_x = np.array([approach.opik.x*au2units for approach in ca_list])
    opik_y = np.array([approach.opik.y*au2units for approach in ca_list])
    scaled_x = np.array([approach.scaled.x*au2units for approach in ca_list])
    scaled_y = np.array([approach.scaled.y*au2units for approach in ca_list])
    mtp_x = np.array([approach.mtp.x*au2units for approach in ca_list])
    mtp_y = np.array([approach.mtp.y*au2units for approach in ca_list])
    focus_factor = np.max([approach.gravFocusFactor for approach in ca_list])
    impact_bool = np.array([approach.impact for approach in ca_list])
    impact_any = np.any(impact_bool)
    if len(ca_list) >= 13:
        kizner_ellipse = data_to_ellipse(kizner_x, kizner_y, n_std, plot_offset,
                                            'kizner', print_ellipse_params, units, sigma_points)
        opik_ellipse = data_to_ellipse(opik_x, opik_y, n_std, plot_offset,
                                            'opik', print_ellipse_params, units, sigma_points)
        scaled_ellipse = data_to_ellipse(scaled_x, scaled_y, n_std, plot_offset,
                                            'scaled', print_ellipse_params, units, sigma_points)
        mtp_ellipse = data_to_ellipse(mtp_x, mtp_y, n_std, plot_offset,
                                            'mtp', print_ellipse_params, units, sigma_points)
    else:
        kizner_ellipse = None
        opik_ellipse = None
        scaled_ellipse = None
        mtp_ellipse = None
    fig, axes = plt.subplots(2, 2, figsize=(9, 9), dpi=150)
    plot_single_bplane(axes[0,0], kizner_x, kizner_y, kizner_ellipse, 'kizner',
                        focus_factor, show_central_body, plot_offset, scale_coords,
                        central_body_radius, units, equal_axis)
    plot_single_bplane(axes[0,1], opik_x, opik_y, opik_ellipse, 'opik',
                        focus_factor, show_central_body, plot_offset, scale_coords,
                        central_body_radius, units, equal_axis)
    plot_single_bplane(axes[1,0], scaled_x, scaled_y, scaled_ellipse, 'scaled',
                        focus_factor, show_central_body, plot_offset, scale_coords,
                        central_body_radius, units, equal_axis)
    plot_single_bplane(axes[1,1], mtp_x, mtp_y, mtp_ellipse, 'mtp',
                        focus_factor, show_central_body, plot_offset, scale_coords,
                        central_body_radius, units, equal_axis)
    patches, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(patches, labels, loc='upper center', ncol=4,
                bbox_to_anchor=(0.5, 1.023), fontsize=10)
    fig.tight_layout()
    event = "Impact" if impact_any else "Close Approach"
    unit = "UTC" if impact_any else "TDB"
    if sigma_points is None:
        t_mean = np.mean(times)
        t_dev = np.std(times)
        t_map_mean = np.mean(map_times)
    else:
        t_mean, t_var = sigma_points.reconstruct(times)
        t_dev = np.sqrt(t_var[0,0])
        t_map_mean, _ = sigma_points.reconstruct(map_times)
    if impact_any:
        t_mean_str = Time(t_mean, format='mjd', scale='tdb').utc.iso
        t_map_mean = Time(t_map_mean, format='mjd', scale='tdb').utc.iso
    else:
        t_mean_str = Time(t_mean, format='mjd', scale='tdb').tdb.iso
        t_map_mean = Time(t_map_mean, format='mjd', scale='tdb').tdb.iso
    t_std = n_std*t_dev
    *_, t_std_str = days_to_dhms(t_std)
    full_str = fr'{t_mean_str} $\pm$ {t_std_str}'
    fig.suptitle(fr"{event} at {full_str} {unit} ({n_std}$\sigma$)", fontsize=14, y=1.07)
    subtitle = f"B-plane map time: {t_map_mean} {unit}"
    plt.text(x=0.5, y=1.026, s=subtitle, fontsize=11, ha="center", transform=fig.transFigure)
    plt.show()
    return None

def _get_bplane_labels(bplane_type, plot_offset, units):
    """
    Get the labels and titles for the B-plane plots.

    Parameters
    ----------
    bplane_type : str
        type of B-plane
    plot_offset : bool
        True to plot the offset from the mean, False to plot the raw data
    units : str
        units of the plot

    Returns
    -------
    x_label : str
        x-axis label
    y_label : str
        y-axis label
    title : str
        plot title
    """
    if bplane_type == 'kizner':
        if plot_offset:
            x_label = fr"B.R - B.R$_{{nom}}$ [{units}]"
            y_label = fr"B.T - B.T$_{{nom}}$ [{units}]"
        else:
            x_label = fr"B.R [{units}]"
            y_label = fr"B.T [{units}]"
        title = "Kizner B-plane"
    elif bplane_type == 'opik':
        if plot_offset:
            x_label = fr"$\xi$ - $\xi_{{nom}}$ [{units}]"
            y_label = fr"$\zeta$ - $\zeta_{{nom}}$ [{units}]"
        else:
            x_label = fr"$\xi$ [{units}]"
            y_label = fr"$\zeta$ [{units}]"
        title = "Öpik B-plane"
    elif bplane_type == 'scaled':
        if plot_offset:
            x_label = fr"B.R$_{{scaled}}$ - B.R$_{{scaled, nom}}$ [{units}]"
            y_label = fr"B.T$_{{scaled}}$ - B.T$_{{scaled, nom}}$ [{units}]"
        else:
            x_label = fr"B.R$_{{scaled}}$ [{units}]"
            y_label = fr"B.T$_{{scaled}}$ [{units}]"
        title = "Scaled Kizner B-plane"
    elif bplane_type == 'mtp':
        if plot_offset:
            x_label = fr"X - X$_{{nom}}$ [{units}]"
            y_label = fr"Y - Y$_{{nom}}$ [{units}]"
        else:
            x_label = f"X [{units}]"
            y_label = f"Y [{units}]"
        title = "Modified Target Plane"
    return x_label, y_label, title

def plot_single_bplane(axis, x_coord, y_coord, ellipse, bplane_type,
                        focus_factor, show_central_body, plot_offset, scale_coords,
                        central_body_radius, units, equal_axis):
    """
    Plot a single B-plane.

    Parameters
    ----------
    axis : matplotlib axis
        axis to plot on
    x_coord : np.ndarray
        x-coordinates of the points
    y_coord : np.ndarray
        y-coordinates of the points
    ellipse : np.ndarray
        ellipse points
    bplane_type : str
        type of B-plane
    focus_factor : float
        scaling factor of the central body B-plane due to gravitational focusing
    show_central_body : bool
        True to show the central body, False to not show it
    plot_offset : bool
        True to plot the offset from the mean, False to plot the raw data
    scale_coords : bool
        True to scale the coordinates to the body's radius, False to leave in km
    central_body_radius : float
        radius of the central body, in km
    units : str
        units of the plot
    equal_axis : bool
        True to make the x and y axes equal, False to leave them as is

    Returns
    -------
    None : NoneType
        None

    Raises
    ------
    ValueError
        Unknown B-plane type
    """
    if bplane_type not in {'kizner', 'opik', 'scaled', 'mtp'}:
        raise ValueError("Unknown B-plane type")
    scale_factor = central_body_radius if scale_coords else 1.0
    if bplane_type in {'scaled', 'mtp'}:
        focus_factor = 1.0
    x_label, y_label, title = _get_bplane_labels(bplane_type, plot_offset, units)
    msize = 5
    malpha = 0.5
    mspec = 'b.'
    rotation = 30
    axis.plot(x_coord/scale_factor, y_coord/scale_factor, mspec,
                ms=msize, alpha=malpha, label="Close Approaches")
    if ellipse is not None:
        lwidth = 1.0
        lalpha = 0.75
        lspec = 'r-'
        axis.plot(ellipse[0,:]/scale_factor, ellipse[1,:]/scale_factor, lspec,
                    lw=lwidth, alpha=lalpha, label="Uncertainty Ellipse")
    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    axis.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useOffset=False)
    axis.tick_params(axis='x', labelrotation=rotation)
    axis.set_title(title, fontsize=12)
    cclr = 'k'
    calpha = 0.5
    clst = '-'
    if scale_coords:
        circle = plt.Circle((0, 0), 1, label="Central Body Extent",
                            color=cclr, fill=True, alpha=0.55*calpha, linestyle=clst)
        circle2 = plt.Circle((0, 0), focus_factor, label="Impact Cross Section",
                            color=cclr, fill=True, alpha=0.4*calpha, linestyle=clst)
    else:
        circle = plt.Circle((0, 0), central_body_radius, label="Central Body Extent",
                            color=cclr, fill=True, alpha=0.55*calpha, linestyle=clst)
        circle2 = plt.Circle((0, 0), focus_factor*central_body_radius, label="Impact Cross Section",
                            color=cclr, fill=True, alpha=0.4*calpha, linestyle=clst)
    if show_central_body:
        axis.add_patch(circle)
        axis.add_patch(circle2)
    if equal_axis:
        axis.axis('equal')
    return axis
