// only adjacent brackets can be swapped
public class MinSwapsToBalanceBracket {
	
	static int swapCount(String s)
	{	
		int open = 0;
		int close = 0;
		int swaps = 0;
		int imbalance = 0;
		
		for(int i = 0;i< s.length();i++)
		{
			if(s.charAt(i) == '[')
			{
				open++;
				if(imbalance > 0)
				{
					swaps += imbalance;
					imbalance--;
				}
			}
			
			else if(s.charAt(i) == ']')
			{
				close++;
				imbalance = close - open;
			}
		}
		
		return swaps;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String s = "[]][][][][";
		System.out.println(swapCount(s));
		
		String s1 = "[]][][]";
		System.out.println(swapCount(s1));
		

	}

}
