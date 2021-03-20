package arrays;

public class BuySellStock5 {
// buy and sell stocks if two transaction is allowed	
	static int maxProfit(int price[])
	{
		int n = price.length;
		int mpist = 0;    // max profit if sold today
		int leastsf = price[0];      // least so far
		int dpl[] = new int[n];
		
		for(int i=1;i<n;i++)
		{
			if(price[i]< leastsf)
			{
				leastsf = price[i];
			}
			
			mpist = price[i] - leastsf;
			
			dpl[i] = Math.max(dpl[i-1], mpist);
		}
		
		int mpibt =0;
		int maxsf = price[n-1];
		
		int dpr[] = new int[n];
		
		for(int i=n-2;i>=0;i--)
		{
			if(price[i] > maxsf)
			{
				maxsf = price[i];
			}
			mpibt = maxsf - price[i];
			
			dpr[i]  =Math.max(dpr[i+1], mpibt);
		}
		
		int op = 0;
		for(int i=0;i<n;i++)
		{
			op = Math.max(op, dpl[i]+dpr[i]);
		}
		return op;
	}

	
	public static void main(String[] args) {
		
		int price[] = {10, 20, 30};
		System.out.println(maxProfit(price));

}
	}
