import {test,expect} from '@playwright/test';
test('missing genotype review and German controls',async({page})=>{
 await page.goto('/');await expect(page.getByRole('heading',{name:'Study overview',exact:true})).toBeVisible();
 await page.getByRole('button',{name:'miHA discovery',exact:true}).click();await page.getByLabel('Search candidates').fill('PTK2B');await page.getByRole('button',{name:'PTK2B',exact:true}).click();
 await expect(page.getByText('No-call: donor absence is unproven',{exact:true})).toBeVisible();
 await page.getByRole('button',{name:'Switch language'}).click();await expect(page.getByRole('heading',{name:'miHA-Analyse',exact:true})).toBeVisible();
});
test('controlled negative assay downgrades functional evidence',async({page})=>{
 await page.goto('/');await page.getByRole('button',{name:'Validation workbench',exact:true}).click();
 await page.getByRole('combobox',{name:'Assay',exact:true}).selectOption('killing');await page.getByLabel('Outcome',{exact:true}).selectOption('negative');await page.getByLabel('Replicates',{exact:true}).fill('3');await page.getByLabel('Controls',{exact:true}).selectOption('adequate');await page.getByLabel('Source / experiment reference').fill('Synthetic controlled negative QA');
 await page.getByRole('button',{name:'Add evidence record'}).click();await expect(page.getByText('conflicting',{exact:true})).toBeVisible();
 await page.getByRole('button',{name:'Study overview',exact:true}).click();await expect(page.getByRole('row').filter({hasText:'MIHA-01'})).toContainText('Binding evidence');
});
test('report exports reproducible inputs and shortlist',async({page})=>{
 await page.goto('/');await page.getByRole('button',{name:'Validation workbench',exact:true}).click();await page.getByRole('button',{name:'Save for review'}).click();
 const download=page.waitForEvent('download');await page.getByRole('button',{name:'Export report'}).click();const d=await download;const stream=await d.createReadStream();const parts:Buffer[]=[];for await(const part of stream!)parts.push(part);const report=JSON.parse(Buffer.concat(parts).toString());expect(report.input_sha256).toMatch(/^[a-f0-9]{64}$/);expect(report.shortlist).toHaveLength(1);expect(report.input.project.synthetic).toBe(true);
});
test('mobile layout and study pair isolation',async({page})=>{
 await page.setViewportSize({width:390,height:844});await page.goto('/');await page.getByLabel('Study pair').selectOption('TX-002');await page.getByRole('button',{name:'miHA discovery',exact:true}).click();await expect(page.getByRole('button',{name:'ARHGAP45',exact:true})).toBeVisible();await expect(page.getByText('Recipient germline provenance is unresolved',{exact:true})).toBeVisible();
 expect(await page.evaluate(()=>document.documentElement.scrollWidth<=window.innerWidth+2)).toBe(true);
});
