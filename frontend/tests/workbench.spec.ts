import {test,expect} from '@playwright/test';

test('explore, inspect evidence and export report',async({page})=>{
 await page.goto('/');
 await expect(page.getByRole('button',{name:'EGFR',exact:true})).toBeVisible();
 await page.getByRole('button',{name:'EGFR',exact:true}).click();
 await expect(page.getByRole('heading',{name:'Why this target appears'})).toBeVisible();
 await expect(page.getByText('Amplification CN approximately 6;', {exact:false})).toBeVisible();
 await page.getByRole('button',{name:'Pathway story'}).click();
 await expect(page.getByRole('heading',{name:'EGFR: observation'})).toBeVisible();
 const downloading=page.waitForEvent('download');
 await page.getByRole('button',{name:'Export report'}).click();
 expect((await downloading).suggestedFilename()).toBe('causalomics-report.json');
 await page.screenshot({path:'test-results/desktop.png',fullPage:true});
});

test('upload CellOmics CSV changes actual analysis',async({page})=>{
 await page.goto('/');
 await page.getByRole('button',{name:'Data workspace'}).click();
 await page.getByLabel('Upload evidence').setInputFiles({name:'example.csv',mimeType:'text/csv',buffer:Buffer.from('gene,tumor_tpm,control_tpm\nBRCA1,10,100\n')});
 await expect(page.getByRole('button',{name:'BRCA1',exact:true})).toBeVisible();
 await expect(page.getByRole('button',{name:'EGFR',exact:true})).toHaveCount(0);
 await expect(page.getByText('1 evidence records')).toBeVisible();
});

test('bad upload preserves dataset and shows error',async({page})=>{
 await page.goto('/');await page.getByRole('button',{name:'Data workspace'}).click();
 await page.getByLabel('Upload evidence').setInputFiles({name:'bad.csv',mimeType:'text/csv',buffer:Buffer.from('gene,tumor_tpm,control_tpm\nEGFR,nan,100\n')});
 await expect(page.getByRole('alert')).toContainText('finite');
 await page.getByRole('button',{name:'Discovery',exact:false}).click();
 await expect(page.getByRole('button',{name:'EGFR',exact:true})).toBeVisible();
});

test('mobile navigation and layout',async({page})=>{
 await page.setViewportSize({width:390,height:844});await page.goto('/');
 await expect(page.getByRole('button',{name:'EGFR',exact:true})).toBeVisible();
 expect(await page.evaluate(()=>document.documentElement.scrollWidth)).toBeLessThanOrEqual(390);
 await page.screenshot({path:'test-results/mobile.png',fullPage:true});
 await page.getByRole('button',{name:'Methods'}).click();
 await expect(page.getByRole('heading',{name:'Evidence weights'})).toBeVisible();
});

test('minimal JSON imports with server defaults',async({page})=>{
 await page.goto('/');await page.getByRole('button',{name:'Data workspace'}).click();
 const data={name:'minimal',evidence:[{id:'1',gene:'TP53',layer:'genomics',disease:'D',tissue:'T',study:'test',source:'test',observation:'fixture',strength:0.5}]};
 await page.getByLabel('Upload evidence').setInputFiles({name:'minimal.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(data))});
 await expect(page.getByRole('button',{name:'TP53',exact:true})).toBeVisible();
 await page.getByRole('button',{name:'Data workspace'}).click();
 await expect(page.getByRole('heading',{name:'Provenance & preparation'})).toBeVisible();
});
