
public class SmallestWindowContainingChOfString {
	
	static void window(String s)
	{
		int n = s.length();
		int count[] = new int[256];
		
		for(int i =0 ;i<n;i++)
		{
			count[(int)s.charAt(i)]++;
		}
		
		int i=0;
		int j = n-1;
		
		while(i<=j) 
		{	
		if(count[s.charAt(i)] > 1)
		{
			count[s.charAt(i)]--;
			i++;
			
		}
		else if(count[s.charAt(j)] > 1)
		{
			count[s.charAt(j)]--;
			j--;
		}
		else
			break;
		}

		System.out.println(s.substring(i, j+1));
		
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String s= "aaab";
		String s1= "aabcbcdbca";
		String s2 = "aaaabbbbbd";
		window(s2);

	}

}
