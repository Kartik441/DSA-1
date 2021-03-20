// Split binary String into substring with equal number of 0s and 1s
public class SplitBinaryStringWithEqual0s1s {
	
static int solve(String str)
{
	int n = str.length();
	int count = 0;
	int ans = 0;
	for(int i=0;i<n;i++)
	{
		if(str.charAt(i) == '0')
			count++;
		else
			count--;
		
		if(count == 0)
			ans++;
	
	}
	return ans;
}
	public static void main(String[] args) {
		
		String str = "0100110101";
		System.out.println(solve(str));
		String str2 = "0111100010";
		System.out.println(solve(str2));

	}

}
