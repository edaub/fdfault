#ifndef COORDCLASSHEADERDEF
#define COORDCLASSHEADERDEF

class coord
{
public:
	coord();
    int get_nx(const int index) const;
    int get_xm(const int index) const;
	int get_xp(const int index) const;
    int get_nx_loc(const int index) const;
	int get_nx_tot(const int index) const;
    int get_xm_loc(const int index) const;
	int get_xp_loc(const int index) const;
    int get_xm_ghost(const int index) const;
    int get_xp_ghost(const int index) const;
	int get_min_loc(const int index) const;
	int get_max_loc(const int index) const;
    void set_nx(const int index, const int new_nx);
    void set_xm(const int index, const int new_xm);
    void set_nx_loc(const int index, const int new_nx);
    void set_xm_loc(const int index, const int new_xm);
    void set_xm_ghost(const int index, const int new_xm);
    void set_xp_ghost(const int index, const int new_xp);
private:
	int nx[3];
    int xm[3];
    int nx_loc[3];
    int xm_loc[3];
    int xm_ghost[3];
    int xp_ghost[3];
};

#endif