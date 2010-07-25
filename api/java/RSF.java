package rsf;

public class RSF{

    public RSF(String[] args){      
        if (args.length == 0){
            args = new String[]{"junk"};
        }  else if (args.length == 1){
            args = new String[]{"junk",args[0]};
        }
        m8r.sf_init(args.length,args);
    }

    public boolean getBool(String key, boolean fallback){
        boolean found = m8r.sf_getbool(key);
        return found;
    }

    public int getInt(String key, int fallback){
        int[] value = new int[]{0};
        boolean found = m8r.sf_getint(key,value);
        if (found) return value[0];
        else return fallback;
    }

    public float getFloat(String key, float fallback){
        float[] value = new float[]{0.0f};
        boolean found = m8r.sf_getfloat(key,value);
        if (found) return value[0];
        else return fallback;
    }

    public String getString(String key, String fallback){
        String value = m8r.sf_getstring(key);
        if (value == null) return fallback;
        else return value;
    }
}
