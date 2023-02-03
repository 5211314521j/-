/*
import java.util.Date;
import java.util.Scanner;

public class Main {
    public static void main(String[] args) {
        //轧制压力模型
        //B:扎件的宽度
        double B;
        //R1:扎件的压扁半径
        double R1;
        //Qp 应力状态系数
        double Qp;
        //K 平均累积平面变形抗力
        double K;
        //Nt 张力影响系数
        double Nt;
        //R
        double R;
        //tb :后张应力
        //tf:前张应力
        double tb;
        double tf;
        //Vo:轧辊线速度
        double Vo;
        double a=3.33;
        double f;
        System.out.println(" 轧制压力模型");
       Scanner sc =new Scanner(System.in);
        System.out.println("请输入后张应力");
         tb = sc.nextDouble();
        System.out.println("请输入前张应力");
        tf =sc.nextDouble();
        System.out.println("请输入平均累积平面变形抗力");
        K=sc.nextDouble();
        Nt = CalculateTheTensionInfluenceCoefficient(a, tb, tf, K);
        System.out.println("张力影响系数为"+Nt);
        //输入轧辊线速度，计算f
        System.out.println("请输入轧辊线速度Vo");
        Vo=sc.nextDouble();
         f = CalculateF(Vo);

        System.out.println("请输入B");
         B =sc.nextDouble();
        System.out.println("请输入H");
        double H =sc.nextDouble();
        System.out.println("请输入h");
        double h =sc.nextDouble();
        System.out.println("请输入Co");
        double Co =sc.nextDouble();
        System.out.println("请输入R");
        R =sc.nextDouble();
        R1 = CalculateR1(B, H, h, Co, R);
         Qp = CalculateQp(R1, H, f);
        System.out.println("变形抗力模型");
        System.out.println("请输入Co");
        double eh=sc.nextDouble();
        System.out.println("请输入R");
        double eH =sc.nextDouble();
         K = Calculatek(eh, eH);
        double P = CalculateP(B, K, R1, H, h, Qp, Nt);
        System.out.println(P);

    }

    private static double CalculateP(double b, double k, double r1, double H, double h, double qp, double nt) {
        //根据公式计算P
        double P =b*k*Math.sqrt(r1*(H-h))*qp*nt;
      return P;
    }

    private static double Calculatek(double eh, double eH) {
        //计算K
        //a1,a2,a3为回归系数
        double a1=-2.992;
        double a2=3.72;
        double a3=0.2451;
        double e1=0.4*eH+0.6*eh;
        double k =a1*Math.pow((e1+a2),a3);
        return k;
    }

    private static double CalculateQp(double R1, double H,double f) {
       //根据公式求QP应力状态系数
        //E为变形抗力
        //
        double E;
        double Qp;
        Qp =1.08+1.79*f*Math.sqrt((R1/H))-1.02*E;
        return Qp;
    }

    private  static  double CalculateR1(double B,double H,double h,double Co,double R){
        //根据公式求R1:扎件的压扁半径
        double R1;
        double P;
        double sum1 =B*(H-h);
        double sum2 =Co*(P/(sum1));
        R1 =R*(1+sum2);
        return  R1;
    }

    private static double CalculateF(double vo) {
        //根据公式求f  vo为轧辊线速度
        double sum1 =3.28*vo;
        double sum2 =4.9*Math.pow(sum1,-0.038);
        double f=0.001*Math.exp(sum2);
        return f;
    }

    private static double CalculateTheTensionInfluenceCoefficient(double a, double tb, double tf, double k) {
        //Nt 求张力影响系数  tb 后张应力  tf前张应力
        double Nt;
        double  a1=(a-1)*tb+tf;
        double a2=a*k;
        Nt =1-(a1/a2);
        return Nt;
    }
}*/
