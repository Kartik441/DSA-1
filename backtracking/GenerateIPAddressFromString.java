package backtracking;

public class GenerateIPAddressFromString {
	// present number of digits in a present block in ip address
	static int digCount = 0;
	
	static void generate(String s, int i, int count, String ip)
	{
		// if whole string is traversed and count is equal to 3 
		if(s.length() == i && count == 3 && digCount!=0)   // digCount for making sure not an output like this
		{                                                  //  11.35.45.
			System.out.println(ip);
			return;
			
		}
		
		// if string is traversed but count is not equal to three
		if(s.length() <= i)
			return;
		
		// adding the character at i in ip address
		ip = ip +s.charAt(i);
		// increment number of digits in block
		digCount++;
		
		// if number in any block exceeds three than return
		if(digCount > 3)
			return;
		
		int x = digCount;
		
		// if first digit of block equals to zero than return
		if(digCount == 1)
		{
			if(s.charAt(i) == '0')
				return;
		}
		
		// choosing to put . OR putting end to this block
		
		ip = ip + ".";
		
		if(digCount == 3) {
		// converting present block into number
		int value = Integer.parseInt(ip.substring(ip.length() - digCount - 1, ip.length()-1));
	
		//  if block value exceeds 255 than return 
		if(value > 255)
			return;
		}
		// since ended the valid block hence digCount = 0 for new block
		digCount = 0;
		generate(s, i+1, count+1, ip);
		
		// backtrack or removing the dot(.)
		ip = ip.substring(0, ip.length()-1);
		
		// retaining the original digCount from x
		digCount = x;
		
		generate(s, i+1, count, ip);
		
		
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String s = "25525511135";
		String s1 = "11111011111";
		String s2 = "1111";
		String s3 = "25505011535";
		generate(s, 0, 0, "");
		digCount = 0;
		
		generate(s1, 0, 0, "");
		digCount = 0;
		
		generate(s2, 0, 0, "");
		digCount = 0;
		
		generate(s3, 0, 0, "");
		digCount = 0;

	}

}
