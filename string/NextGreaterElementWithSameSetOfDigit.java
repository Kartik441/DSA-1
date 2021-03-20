

public class NextGreaterElementWithSameSetOfDigit {
	
	static void next(String s)
	{
		int n = s.length();
		int i = n-2;
		for(i = n-2;i>=0;i--)
		{
			if(s.charAt(i) < s.charAt(i+1))
			{
				break;
			}
				
		}
		
		if(i == -1)
		{
			
			System.out.println("Not Possible");
		}
		else
		{
			s = swap(s, i, n-1);
			// now sort the substring
			// reversing the substring will give sorted string
			String sorted = reverse(s.substring(i+1));
			
			s = s.substring(0, i+1) + sorted;
			
			System.out.println(s);
		}
		
	}
	
	static String swap(String s, int i, int j)
	{
		char a[] = s.toCharArray();
		char temp = a[i];
		a[i] = a[j];
		a[j] = temp;
		
		return String.valueOf(a);
		
		
	}
	
	static String reverse(String s)
	{
		int n=s.length();
		char arr[] = s.toCharArray();
		int i=0;                                                           //  0;
		int j=n-1;
		while(i<j)
		{
			char temp = arr[i];
			arr[i] = arr[j];
			arr[j] = temp;
			i++;
			j--;
		}
		
		return String.valueOf(arr);
	}

	public static void main(String[] args) {
		
		String s ="534976";
		String s1 = "218765";
		String s2 = "1234";
		String s3 = "4321";
		next(s);
		next(s1);
		next(s2);
		next(s3);

	}

}
