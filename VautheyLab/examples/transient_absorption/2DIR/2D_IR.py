from VautheyLab.transient_absorption import twoDIR

t = twoDIR('data.p2dat',
           delay=0.5,               # population delay T2 in ps
           xlim=[2170, 2270],       # x-axis limit
           ylim=[2170, 2270],       # y-axis limit
           flip=False,              # flip exc. / probe axis
           lines=True,              # whether to plot black lines on contours
           nlevels=50,              # levels
           scale=[-0.03, 0.03])     # colorbar scale

#t.get_anharmonicity(range_probe=[2235, 2260], range_pump=[2180, 2300])
t.show()