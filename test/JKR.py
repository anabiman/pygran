font = {'family' : 'normal',
                          
        'size'   : 22}

rc('font',**font)

E = 7.1e7; R = 1e-5; gamma = 20.0;

a = arange(0, R*0.1, R*0.001)

delta = a**2 / R - sqrt(2 * pi * gamma * a / E)

F = 4.0 / 3.0 * E / R * a**3.0 - sqrt(8.0 * pi * gamma * E * a**3.0)

ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
