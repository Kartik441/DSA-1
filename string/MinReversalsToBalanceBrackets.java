
public class MinReversalsToBalanceBrackets {
	
	//can be solved recursively O(2^n)

	// The idea is to first remove all balanced part of expression
	// O(n) solution
	// just count number of unbalanced open and close brackets
	// and add ceil value of them after dividing with 2
	static int min(String s)
	{
		int n = s.length();
		int open = 0;       // current open that are unbalanced
		int close = 0;      // current close that are unbalanced
		
		if(s.length() % 2 !=0)
			return -1;
		
		for(int i=0;i<n;i++)
		{
			if(s.charAt(i) == '{')
			{
				open++;
			}
			else
			{
				if(open == 0)
					close++;
				else
					open--;
			}
		}
		
		// after this loop open and close will contain unbalanced brackets
//		int ans = (close/2) + (open/2);
//		
//		// for case: }{
//		
//		if(close % 2 != 0)     // if close is odd
//			ans += 2;
		
		int ans = 0;
		ans += Math.ceil((double)close/2);
		ans += Math.ceil((double)open/2);
		
		
		return ans;
		
		
		
	}

	public static void main(String[] args) {
		
		String s = "}{{}}{{{";
		String s1 = "}{}}";
		String s2 = "}}{}{{";     // ??????
		//System.out.println(Math.ceil(0.5));
		System.out.println(min(s));
		System.out.println(min(s1));
		System.out.println(min(s2));
		String s3 = "{{{{";
		System.out.println(min(s3));

	}

}
