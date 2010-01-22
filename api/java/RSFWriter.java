package rsf;

import java.io.*;
import java.nio.ByteOrder;
import edu.mines.jtk.io.ArrayOutputStream;

public class RSFWriter{

    public static void writeHeader(RSFHeader header,String fullPath){

            try{
                    PrintWriter pw = new PrintWriter(new FileOutputStream(header.name));
                    String path = System.getenv("DATAPATH");
                    pw.println("java programs warning, this file may not be reproducibile");
                    pw.println();
                    pw.println("in=\""+fullPath+"\"");
                    pw.println("data_format="+"\""+header.format+"\"");
                    pw.println("esize="+header.esize);
                    for(int i = 0, j=1; i < header.ndims; ++i, ++j){
                        pw.println("n"+j+"="+header.n[i]);
                        pw.println("d"+j+"="+header.deltas[i]);
                        pw.println("o"+j+"="+header.origins[i]);
                        pw.println("label"+j+"=\""+header.labels[i]+"\"");
                        pw.println("unit"+j+"=\""+header.units[i]+"\"");
                    }
                    pw.close();
            } catch (Exception e){

                    System.out.println(e);
            }
                 
    }

    public static void writeRSF(RSFHeader header, float[] data){
            String fullPath = System.getenv("DATAPATH")+ header.name+"@";
            writeHeader(header,fullPath);
            try{
                    ArrayOutputStream aos = new ArrayOutputStream(fullPath,ByteOrder.nativeOrder());

                    aos.writeFloats(data);

                    aos.close();
            } catch (Exception e){
                    System.out.println(e);
                    }

    }

    public static void writeRSF(RSFHeader header, float[][] data){
            String fullPath = System.getenv("DATAPATH")+ header.name+"@";
            writeHeader(header,fullPath);
            try{
                   ArrayOutputStream aos = new ArrayOutputStream(fullPath,ByteOrder.nativeOrder());

                    aos.writeFloats(data);

                    aos.close();
            } catch (Exception e){
                    System.out.println(e);
            }
    }
    public static void writeRSF(RSFHeader header, float[][][] data){
            String fullPath = System.getenv("DATAPATH")+ header.name+"@";
            writeHeader(header,fullPath);
            try{
                   ArrayOutputStream aos = new ArrayOutputStream(fullPath,ByteOrder.nativeOrder());

                    aos.writeFloats(data);

                    aos.close();
            } catch (Exception e){
                    System.out.println(e);
            }
    }


}

