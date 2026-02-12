from VautheyLab.transient_absorption import Compare_Overviews

o = Compare_Overviews(files = ['', ''],
                      scatter=[(), ()],
                      MA=[True, True],
                      MA_npoints=[10, 10])
o.show()