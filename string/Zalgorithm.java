
public class Zalgorithm {

	static int[] zFunction(String s)
	{
		int z[] = new int[s.length()];
		// at i index storage is
		// the length of the longest substring starting from i index which is also a proper prefix
		
		int l = 0;
		int r = 0;
		for(int i=1;i<s.length();i++)
		{
			if(i <= r)
			{
				// give a kickstart
				z[i] = Math.min(r - i + 1, z[i-l]);
			}
			
			while(s.charAt(z[i]) == s.charAt(i + z[i]))
			{
				z[i]++;
			}
			
			if(i + z[i] - 1 > r)
			{
				// update the info
				l = i;
				r = i + z[i] - 1;
			}
		}
		for(int i:z)
			System.out.print(i+" ");
		return z;
	}
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		zFunction("ababadababae");
	}

}
