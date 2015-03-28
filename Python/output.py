from __future__ import division, print_function

class output(object):
    def __init__(self, name, field, tm = 0, tp = 0, ts = 1, xm = 0, xp = 0, xs = 1, ym = 0,
                 yp = 0, ys = 1, zm = 0, zp = 0, zs = 1):
        "Initializes ouput unit"
        assert type(name) is str, "output name must be a string"
        assert (field == "vx" or field == "vy" or field == "vz" or field == "sxx" or field == "sxy"
                or field == "sxz" or field == "syy" or field == "syz" or field == "szz" or field == "Ux"
                or field == "Uy" or field == "Vx" or field == "Vy" or field == "U" or field == "V"), "Incorrect field name"
        assert (tm >= 0 and tp >= 0 and tp >= tm and ts > 0), "bad time limits"
        assert (xm >= 0 and xp >= 0 and xp >= xm and xs > 0), "bad x limits"
        assert (ym >= 0 and yp >= 0 and yp >= ym and ys > 0), "bad y limits"
        assert (zm >= 0 and zp >= 0 and zp >= zm and zs > 0), "bad z limits"
        
        self.name = name
        self.field = field
        self.tm = tm
        self.tp = tp
        self.ts = ts
        self.xm = xm
        self.ym = ym
        self.zm = zm
        self.xp = xp
        self.yp = yp
        self.zp = zp
        self.xs = xs
        self.ys = ys
        self.zs = zs

    def get_name(self):
        "Returns name"
        return self.name

    def set_name(self, name):
        "Set output unit name"
        assert type(name) is str, "output name must be a string"
        self.name = name

    def get_field(self):
        "Returns output field"
        return self.field

    def set_field(self, field):
        assert (field == "vx" or field == "vy" or field == "vz" or field == "sxx" or field == "sxy"
                or field == "sxz" or field == "syy" or field == "syz" or field == "szz" or field == "Ux"
                or field == "Uy" or field == "Vx" or field == "Vy" or field == "U" or field == "V"), "Incorrect field name"
        self.field = field

    def get_tm(self):
        "Returns minimum time step for output"
        return self.tm

    def get_tp(self):
        "Returns maximum time step for output"
        return self.tp

    def get_ts(self):
        "Returns time stride for output"
        return self.ts

    def get_time_indices(self):
        "Returns all index info for time output as (tm, tp, ts)"
        return (self.tm, self.tp, self.ts)

    def get_xm(self):
        "Returns minimum x grid point for output"
        return self.xm

    def get_xp(self):
        "Returns maximum x grid point for output"
        return self.xp

    def get_xs(self):
        "Returns x stride for output"
        return self.xs

    def get_x_indices(self):
        "Returns all index info for x output as (xm, xp, xs)"
        return (self.xm, selt.xp, self.xs)

    def get_ym(self):
        "Returns minimum y grid point for output"
        return self.ym

    def get_yp(self):
        "Returns maximum y grid point for output"
        return self.yp

    def get_ys(self):
        "Returns y stride for output"
        return self.ys

    def get_y_indices(self):
        "Returns all index info for y output as (ym, yp, ys)"
        return (self.ym, selt.yp, self.ys)

    def get_zm(self):
        "Returns minimum z grid point for output"
        return self.zm

    def get_zp(self):
        "Returns maximum z grid point for output"
        return self.zp

    def get_zs(self):
        "Returns z stride for output"
        return self.zs

    def get_z_indices(self):
        "Returns all index info for z output as (zm, zp, zs)"
        return (self.zm, selt.zp, self.zs)

    def set_tm(self, tm):
        "Sets minimum time index for output"
        assert tm >=0 and tm <= tp, "tm must be positive and not greater than tp"
        self.tm = int(tm)

    def set_tp(self, tp):
        "Sets maximum time index for output"
        assert tp >=0 and tm <= tp, "tp must be positive and not less than tm"
        self.tp = int(tp)

    def set_ts(self, ts):
        "Sets t stride for output"
        assert ts > 0, "ts must be positive"
        self.ts = int(ts)

    def set_time_indices(self, tm, tp, ts):
        "Sets all time indices"
        assert (tm >= 0 and tp >= 0 and tp >= tm and ts > 0), "bad time limits"
        self.set_tm(tm)
        self.set_tp(tp)
        self.set_ts(ts)

    def set_xm(self, xm):
        "Sets minimum x index for output"
        assert xm >=0 and xm <= xp, "xm must be positive and not greater than xp"
        self.xm = int(xm)

    def set_xp(self, xp):
        "Sets maximum x index for output"
        assert xp >=0 and xm <= xp, "xp must be positive and not less than xm"
        self.xp = int(xp)

    def set_xs(self, xs):
        "Sets x stride for output"
        assert xs > 0, "xs must be positive"
        self.xs = int(xs)

    def set_x_indices(self, xm, xp, xs):
        "Sets all x indices"
        assert (xm >= 0 and xp >= 0 and xp >= xm and xs > 0), "bad x limits"
        self.set_xm(xm)
        self.set_xp(xp)
        self.set_xs(xs)

    def set_ym(self, ym):
        "Sets minimum y index for output"
        assert ym >=0 and ym <= yp, "ym must be positive and not greater than yp"
        self.ym = int(ym)

    def set_yp(self, yp):
        "Sets maximum y index for output"
        assert yp >=0 and ym <= yp, "yp must be positive and not less than ym"
        self.yp = int(yp)

    def set_ys(self, ys):
        "Sets y stride for output"
        assert ys > 0, "ys must be positive"
        self.ys = int(ys)

    def set_y_indices(self, ym, yp, ys):
        "Sets all y indices"
        assert (ym >= 0 and yp >= 0 and yp >= xm and ys > 0), "bad y limits"
        self.set_ym(ym)
        self.set_yp(yp)
        self.set_ys(ys)

    def set_zm(self, zm):
        "Sets minimum z index for output"
        assert zm >=0 and zm <= zp, "zm must be positive and not greater than zp"
        self.zm = int(zm)

    def set_zp(self, zp):
        "Sets maximum z index for output"
        assert zp >=0 and zm <= zp, "zp must be positive and not less than zm"
        self.zp = int(zp)

    def set_zs(self, zs):
        "Sets z stride for output"
        assert zs > 0, "zs must be positive"
        self.zs = int(zs)

    def set_z_indices(self, zm, zp, zs):
        "Sets all z indices"
        assert (zm >= 0 and zp >= 0 and zp >= zm and zs > 0), "bad z limits"
        self.set_zm(zm)
        self.set_zp(zp)
        self.set_zs(zs)

    def write_input(self,f):
        "Writes output unit to file"
        f.write(self.name+"\n")
        f.write(self.field+"\n")
        f.write(str(self.tm)+" "+str(self.tp)+" "+str(self.ts)+"\n")
        f.write(str(self.xm)+" "+str(self.xp)+" "+str(self.xs)+"\n")
        f.write(str(self.ym)+" "+str(self.yp)+" "+str(self.ys)+"\n")
        f.write(str(self.zm)+" "+str(self.zp)+" "+str(self.zs)+"\n")
        
    def __str__(self):
        return ("Output unit '"+self.name+"': field = "+self.field+", tm = "+str(self.tm)+", tp = "+str(self.tp)+
                ", ts = "+str(self.ts)+", xm = "+str(self.xm)+", xp = "+str(self.xp)+
                ", xs = "+str(self.xs)+"\nym = "+str(self.ym)+", yp = "+str(self.yp)+
                ", ys = "+str(self.ys)+", zm = "+str(self.zm)+", zp = "+str(self.zp)+
                ", zs = "+str(self.zs))
