
class NekoVertex
{
    public:

    int id;
    double x;
    double y;
    double z;
   

    NekoVertex(int idv, double xv, double yv, double zv)
    {
        id = idv; x = xv; y = yv; z = zv;
    }

    void print() const
    {
        Info << id << " " << x << " " << y << " " << z << nl;
    }
};

class NekoHex
{
    public:

    int id;
    std::vector<NekoVertex> vertices;


    NekoHex(int idv, std::vector<NekoVertex> v)
    {
       id = idv; vertices = v;

    }

    void write(std::ofstream & fs) const
    {
        fs.write(reinterpret_cast<const char*>(&id), sizeof id);
        
        for(auto & vi : vertices)
        {
            fs.write(reinterpret_cast<const char*>(&vi.id), sizeof vi.id);
            fs.write(reinterpret_cast<const char*>(&vi.x), sizeof vi.x);
            fs.write(reinterpret_cast<const char*>(&vi.y), sizeof vi.y);
            fs.write(reinterpret_cast<const char*>(&vi.z), sizeof vi.z);
        }
    }

    // faces in terms of vertex indices
    // 1-based!!
    std::vector<int> face(int i) const
    {
        if (i == 1)
        {
            return std::vector<int>({1, 5, 8, 4});
        }
        else if (i == 2)
        {
            return std::vector<int>({2, 6, 7, 3});
        }
        else if (i == 3)
        {
            return std::vector<int>({1, 2, 6, 5});
        }
        else if (i == 4)
        {
            return std::vector<int>({4, 3, 7, 8});
        }
        else if (i == 5)
        {
            return std::vector<int>({1, 2, 3, 4});
        }
        else if (i == 6)
        {
            return std::vector<int>({5, 6, 7, 8});
        }
        else
        {
            FatalErrorIn("NekoHex::face()")
               << abort(FatalError);

        }

    }
};

class NekoZone
{
    public:

    int e;
    int f;
    int p_e;
    int p_f;
    std::array<int, 4> glb_pt_ids;
    int type;

    NekoZone(int ev, int fv, int p_ev, int p_fv,
            std::array<int, 4> idsv, int typev)
    {
       e = ev; f = fv; p_e = p_ev; p_f = p_fv;
       glb_pt_ids = idsv;
       type = typev;

    }

    void write(std::ofstream & fs)
    {
        fs.write(reinterpret_cast<const char*>(&e), sizeof e);
        fs.write(reinterpret_cast<const char*>(&f), sizeof f);
        fs.write(reinterpret_cast<const char*>(&p_e), sizeof p_e);
        fs.write(reinterpret_cast<const char*>(&p_f), sizeof p_f);
        
        for(auto & pi : glb_pt_ids)
        {
            fs.write(reinterpret_cast<const char*>(&pi), sizeof pi);
        }
        fs.write(reinterpret_cast<const char*>(&type), sizeof type);
    

    }
};

