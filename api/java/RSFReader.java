package rsf;

import java.io.*;
import java.util.Scanner;
import java.util.ArrayList;
import java.util.regex.*;

import java.nio.ByteOrder;
import edu.mines.jtk.io.ArrayInputStream;

public class RSFReader{

        /* This is a demo of how to use the interface 
    public static void main(String[] args){

            RSFHeader header = readHeader(args[0]);

            System.out.printf("Header for %s is: \n", args[0]);
            System.out.println(header);

            float[][] data = readBinary2D(header);
            System.out.println("My data is: ");
            for(int i = 0; i < data.length; ++i){
                printFloats(data[i]);
            }
            header.setName("javatest.rsf");

            RSFWriter.writeRSF(header,data);

    }
    */

    public static void printFloats(float[] f){
        for(int i = 0; i < f.length; ++i){
            System.out.printf("%g ",f[i]);
        }
        System.out.print("\n");
    }

    public static float[] readBinary(RSFHeader header){
            try{
                    int n = 1;
                    /* Get total number of elements */
                    for(int i = 0; i < header.ndims; ++i){
                        n *= header.n[i];
                    }
                    float[] data = new float[n]; 
                    ArrayInputStream ais = new ArrayInputStream(header.path, ByteOrder.nativeOrder());
                    ais.readFloats(data);
                    ais.close();
                    return data;
            } catch (Exception e){
                    System.out.println(e);
                    return null;
            }
    }
    public static float[][] readBinary2D(RSFHeader header){
            try{
                    if (header.ndims != 2){
                        throw new Exception("Bad number of dimensions");
                    }
                    int n1 = header.n[0];
                    int n2 = header.n[1];
                    float[][] data = new float[n2][n1]; 
                    ArrayInputStream ais = new ArrayInputStream(header.path, ByteOrder.nativeOrder());
                    ais.readFloats(data);
                    ais.close();
                    return data;
            } catch (Exception e){
                    System.out.println(e);
                    return null;
            }
    }
    public static float[][][] readBinary3D(RSFHeader header){
            try{
                    if (header.ndims != 3){
                        throw new Exception("Bad number of dimensions");
                    }
                    int n1 = header.n[0];
                    int n2 = header.n[1];
                    int n3 = header.n[2];
                    float[][][] data = new float[n3][n2][n1]; 
                    ArrayInputStream ais = new ArrayInputStream(header.path, ByteOrder.nativeOrder());
                    ais.readFloats(data);
                    ais.close();
                    return data;
            } catch (Exception e){
                    System.out.println(e);
                    return null;
            }
    }
 
    public static RSFHeader readHeader(String filename){
            RSFHeader header = new RSFHeader();
            header.setName(filename);
            try{
                    Scanner s = null;
                    if (filename.equals("in")){
                        s = new Scanner(new InputStreamReader(System.in));
                    } else {
                        s = new Scanner(new FileInputStream(filename));
                    }
                    while(s.hasNextLine()){
                        String line = s.nextLine();
                        if( line.contains("=")){
                                line = line.trim();
                                line = line.replace("\t","");
                                line = line.replace("\n","");
                                String[] split = line.split("=");
                                String key = split[0]; String val = split[1];
                                try{
                                        if(key.startsWith("n")){
                                            String temp = key.substring(1);
                                            Integer dim = Integer.valueOf(temp);
                                            Integer elements = Integer.valueOf(val);
                                            header.setN(dim,elements);
                                        } else if (key.startsWith("data_format")){
                                            header.setFormat(val.replace("\"",""));
                                        } else if (key.startsWith("in")){
                                            header.setPath(val.replace("\"",""));
                                        } else if (key.startsWith("d")){
                                            String temp = key.substring(1);
                                            Integer dim = Integer.valueOf(temp);
                                            Float delta = Float.valueOf(val);
                                            header.setDelta(dim,delta);
                                        } else if (key.startsWith("o")){
                                            String temp = key.substring(1);
                                            Integer dim = Integer.valueOf(temp);
                                            Float origin = Float.valueOf(val);
                                            header.setOrigin(dim,origin);
                                        } else if (key.startsWith("label")){
                                            String temp = key.substring(5);
                                            Integer dim = Integer.valueOf(temp);
                                            header.setLabel(dim,val.replace("\"",""));
                                        } else if (key.startsWith("unit")){
                                            String temp = key.substring(4);
                                            Integer dim = Integer.valueOf(temp);
                                            header.setUnit(dim,val.replace("\"",""));
                                        }
                                } catch (ArrayIndexOutOfBoundsException e){
                                        System.out.printf("WARNING: CAUGHT ARRAYOUTOFBOUNDS! TRYING TO ADD %S WITH %S VALUE\n",key,val);                                
                                        System.out.println("This may be ok, if the dimension is greater than 3!");
                                }
                         }
                    }
                    s.close();
            } catch (Exception e){
                    System.out.println(e);
            }
          return header;
    }
}
