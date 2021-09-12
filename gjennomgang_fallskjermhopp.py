from math import exp

from vpython import *


def calc_rho(z: float) -> float:
    """
    Funksjon som regner ut lufttettheten for en gitt høyde.
    Finner tallene på: https://www.grc.nasa.gov/WWW/k-12/airplane/atmosmet.html

    Args:
        z: float - Høyde i meter over havet.
    Returns:
        rho: float - lufttetthet
    """
    if z > 25000:
        T = -131.21  + 0.00299 * z
        p = 2.488 * ((T + 273.1) / 216.6) ** -11.388
    elif z <= 25000 and z >= 11000:
        T = -56.46
        p = 22.65 * exp(1.73 - 0.000157 * z)
    else:
        T = 15.04 - .00649 * z
        p = 101.29 * ((T + 273.1) / 288.08) ** 5.265
    return p / (.2869 * (T + 273.1))

# 3D Animasjon
sky = box(pos=vec(0, -100, 2e4), size=vec(3000, 0.1, 4e4), color=color.blue)
felix = sphere(vel=vec(0,0,0), radius=100)

# Plot
pos_graf = graph(title="posisjon mot tid", xtitle="t (s)", ytitle="z (m)", fast=False)
pos_curve = gcurve(graph=pos_graf)
vel_graf = graph(title="hastighet mot tid", xtitle="t (s)", ytitle="vz (m/s)", fast=False)
vel_curve = gcurve(graph=vel_graf)
energi_graf = graph(title="Energi mot tid", xtitle="t (s)", ytitle="E (J)", fast=False)
pot_curve = gcurve(graph=energi_graf, color=color.blue, label="pot")
kin_curve = gcurve(graph=energi_graf, color=color.red, label="kin")
mek_curve = gcurve(graph=energi_graf, color=color.green, label="mek")
fnet_graf = graph(title="Netto kraft mot tid", xtitle="t (s)", ytitle="Fnet (N)", fast=False)
fnet_curve = gcurve(graph=fnet_graf)
effekt_graf = graph(title="Luftmotstandens effekt mot tid", xtitle="t (s)", ytitle="P (W)", fast=False)
effekt_curve = gcurve(graph=effekt_graf)

# initialbetingelser og konstanter
pos = vec(0, 0, 38969)
vel = vec(0, 0, 0)
m = 70
g = 9.81
D = 0.5

# Overflateareal ~ høyde * skulderbredde
A = 0.8325

# tidsbetingelser
t = 0
dt = 0.1
tmax = 700

# itererer opp til 700 sekunder, eller høyden er null.
while t < tmax and pos.z >= 0:

    # Bestemmer animasjonshastigheten 1 / dt ~ real time * 100
    rate(1/dt * 100)

    t += dt

    # Regner ut kreftene
    Fg = - m * g * vec(0, 0, 1)
    Fr = 0.5 * D * calc_rho(pos.z) * A * vel.z ** 2 * vec(0, 0, 1)

    # Summen av kreftene
    Fnet = Fg + Fr

    # Lagrer forrige posisjon
    forrige_pos = pos

    # Euler-cromers metode for å oppdatere fart og posisjon
    vel = vel + (Fnet / m) * dt
    pos = pos + vel * dt

    # Regner ut potensiell-, kinetisk- og mekanisk energi
    pot = m * g * pos.z
    kin = 0.5 * m * vel.z ** 2
    mek = pot + kin

    # Infinitesimalt tidssteg
    ds = forrige_pos.z - pos.z

    # Effekt = Kraft * strekning / tid
    effekt = Fr.z * ds / dt

    # Oppdaterer 3D animasjonens posisjon og fart
    felix.pos = pos
    felix.vel = vel

    # Oppdaterer plottene for hver iterering
    pos_curve.plot(t, pos.z)
    vel_curve.plot(t, vel.z)
    pot_curve.plot(t, pot)
    kin_curve.plot(t, kin)
    mek_curve.plot(t, mek)
    fnet_curve.plot(t, Fnet.z)
    effekt_curve.plot(t, effekt)
