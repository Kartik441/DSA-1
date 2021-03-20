package arrays;

public class BuySellStock3 {
	
	// buy and sell stock if infinite transactions are allowed but with a transaction fee on every sell	

	static int maxProfit(int []price, int fee)
	{
		
		int obsp = -price[0];    // old bought state profit
		int ossp = 0;            // old sold state profit
		
		for(int i=1; i<price.length;i++)
		{
			int nbsp = 0;
			int nssp = 0;
			
			// for making bsp
			if(ossp - price[i] > obsp)   // if buy a stock today
			{
				nbsp = ossp - price[i];
			}
			else
			{
				nbsp = obsp;
			}
			
			// for making ssp
			if(obsp + price[i] - fee > ossp)
			{
				nssp = obsp + price[i] - fee;
			}
			else
			{
				nssp =ossp;
			}
			
			obsp = nbsp;
			ossp = nssp;

		}
		return ossp;
	}

	public static void main(String[] args) {
		
		int price[] = {10, 20, 30};
		System.out.println(maxProfit(price, 2));

	}

}
