import java.util.HashMap;

// integer to roman and vice versa from 1 to 3999
public class IntegerToRoman {
	
	static String ItoR(int x)
	{
		String units[] = {"", "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX"};
		String tens[] = {"", "X", "XX", "XXX", "XL", "L", "LX", "LXX", "LXXX", "XC"};
		String hundreds[] = {"", "C", "CC", "CCC", "CD", "D", "DC", "DCC", "DCCC", "CM"};
		String thousands[] = {"", "M", "MM", "MMM"};
		
		return thousands[x/1000] + hundreds[(x%1000)/100] + tens[(x%100) / 10] + units[x%10];
		
	}
	
	static int RtoI(String s)
	{
		// roman number are additive and subtractive in nature
		// III = 3
		// subtraction happens when we have a bigger value in front a smaller value(from left to right)
		HashMap<Character , Integer> map = new HashMap<>();
		map.put('I', 1);
		map.put('V', 5);
		map.put('X', 10);
		map.put('L', 50);
		map.put('C', 100);
		map.put('D', 500);
		map.put('M', 1000);
		
		int ans = 0;
		for(int i=0;i<s.length();i++)
		{
			if(i>0 && map.get(s.charAt(i)) > map.get(s.charAt(i-1)))
			{
				// to undo the previous addition we subtract once and then we subtract again
				// due to subtractive nature. Hence, we subtract 2*(i-1)
				ans += map.get(s.charAt(i)) - 2 * map.get(s.charAt(i-1));
			}
			else
			{
				ans += map.get(s.charAt(i));
			}
		}
		
		return ans;
	}

	public static void main(String[] args) {
		System.out.println(ItoR(3999));
		System.out.println(ItoR(248));
		System.out.println(ItoR(1904));
		System.out.println(ItoR(40));
		
		System.out.println(RtoI("MCMIV"));
		System.out.println(RtoI("IX"));
		System.out.println(RtoI("XL"));
		System.out.println(RtoI("CCXLVIII"));

	}

}
