package arrays;

public class BuySellStock4 {
	
	// infinite buy and sell allowed but after a cooldown after each sell
	
	static int maxProfit(int price[])
	{	
		int obsp = -price[0];
		int ossp = 0;
		int ocsp = 0;
		
		for(int i=1;i<price.length;i++)
		{
			int nbsp = 0;
			int nssp =0;
			int ncsp = 0;
			
			//for making bsp
			if(ocsp - price[i] > obsp)
			{
				nbsp = ocsp - price[i];
			}
			else
			{
				nbsp = obsp;
			}
			
			// for making ssp
			if(obsp + price[i] > ossp)
			{
				nssp = obsp + price[i];
			}
			else
			{
				nssp = ossp;
			}
			
			// for making cooldown state
			
			if(ossp > ocsp)
			{
				ncsp = ossp;
			}
			else
			{
				ncsp = ocsp;
			}
			obsp = nbsp;
			ossp = nssp;
			ocsp = ncsp;
			
		}
		
		return ossp;
	}

	public static void main(String[] args) {
		
		int price[] = {10, 20, 30};
		System.out.println(maxProfit(price));

	}

}
