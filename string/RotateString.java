
public class RotateString {
	
	static String leftRotation(String s,int d)
	{
		String ans=s.substring(d)+s.substring(0,d);
		return ans;
	}
	
	static String rightRotation(String s,int d)
	{
		
		String ans=s.substring(s.length()-d)+s.substring(0,s.length()-d);
		return ans;
	}
	
	static boolean isRotated(String s1,String s2)
	{
		String s3=s1+s1;
		return s3.contains(s2);
	}

	public static void main(String[] args) {
		String s="HelloHero";
		String s1=leftRotation(s,4);
		String s2=rightRotation(s,4);
		System.out.println("Left rotation: "+leftRotation(s,4));
		System.out.println("Right rotation: "+rightRotation(s,4));
		System.out.println(isRotated("mypencil","encilmyp"));
		

	}

}
