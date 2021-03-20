
public class ATOI {
	// This must be improved for covering all test cases
	static int atoi(String s)
	{
		int n = s.length();
		int i=0;
		int no = 0;
		while(i<n)
		{
			int digit = s.charAt(i) - '0';
			if(digit <= 9)
			{
				no = no*10 + digit;
			}
			else
			{
				System.out.println("NOt possible");
				return Integer.MIN_VALUE;
			}
			i++;
		}
		return no;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String s = "078978340";
		System.out.println(atoi(s));

	}

}
