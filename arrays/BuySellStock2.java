package arrays;

public class BuySellStock2 {
	
	// buy and sell stock if infinite transactions are allowed
	
	static int maxProfit(int price[])
	{
		// idea is to collect the profit that is increasing (count value of all upstrokes)
		int n = price.length;
		
		int sd=0;
		int bd=0;
		int profit = 0;
		
		for(int i=1; i<n;i++)
		{	
			if(price[i] >= price[i-1])
			{  
				sd++;
			}
			else
			{
				// if next day profit is going to decrease than collect the profit till now
				profit += price[sd] - price[bd];
				// both buy and sell will increase as long as value is decreasing
				bd = sd = i;
			}
		}
		// after loop collect one more last stroke
		profit += price[sd] - price[bd];
		return profit;
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////
	//In a stock market, there is a product with its infinite stocks. The stock prices are given for
	// N days,
	//where arr[i] denotes the price of the stock on the ith day. There is a rule that a customer can
	//buy at
	// most i stock on the ith day. If the customer
	//has an amount of k amount of money initially, find out the maximum number of stocks a customer can
	//buy.
	
	// idea is to get maximum stocks which are cheapest. So we can sort the array in increasing 
	// order. But we need to maintain order of days as well while sorting
	
	//The idea is to use greedy approach, where we should start buying product from the day when the 
	// stock price is least and so on.
	// So, we will sort the pair of two values i.e { stock price, day } according to the stock price,
	// and if stock prices are same, then
	// we sort according to the day.

	public static void main(String[] args) {
		
		int price[] = {3, 1, 2, 3};
		System.out.println(maxProfit(price));
		
		int a1[] = {100, 180, 260, 310, 40, 535, 695};
		System.out.println(maxProfit(a1));
		
		int a2[] = {7, 1, 5, 3, 6, 4};
		System.out.println(maxProfit(a2));

	}

}
